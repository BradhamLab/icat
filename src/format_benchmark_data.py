import pandas as pd
import os
import numpy as np
from scanpy import api as sc
from downstream.src.analysis import utils as dutils
from icat.src import utils
import re

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        adatas = []
        genes = set()
        # lines in CelSeq benchmark data
        lines = ['H1975', 'H2228', 'HCC827']
        # cell mixtures are mixtures of 9 cells
        line_to_mixture = {lines[0]:[9,0,0],
                           lines[1]:[0,9,0],
                           lines[2]:[0,0,9]}
        bench_regex = re.compile('^(.*?)\.')
        for x, meta in zip(snakemake.input['counts'], snakemake.input['meta']):
            count = pd.read_csv(x, index_col=0)
            obs = pd.read_csv(meta, index_col=0)
            # flag for dataset
            bench = os.path.basename(x[:bench_regex.search(x).end() - 1])
            obs['benchmark'] = bench
            obs['index'] = obs.apply(lambda x: '-'.join([x.benchmark, x.name]),
                                     axis=1)
            obs.set_index('index', inplace=True)

            if len(genes) == 0:
                genes = set(count.index.values)
            else:
                genes = genes.intersection(count.index.values)
            
            adatas.append(sc.AnnData(X=count.T.values, obs=obs,
                                     var=pd.DataFrame(count.index.values,
                                                      index=count.index.values,
                                                      columns=['gene'])))
        genes = list(genes)
        # subset genes down to shared genes
        for i in range(len(adatas)):
            adata = adatas[i][:, genes]
            # convert number of cells of each cell type in a mixture to a
            # mixture id. Remove mixtures with less than 9 cells to avoid
            # unclear identity.
            if all([x in adata.obs.columns for x in lines]):
                adata.obs['n_cells'] = adata.obs[lines].apply(lambda x: sum(x),
                                                              axis=1)
                adata.obs['mixture'] = adata.obs[lines].apply(lambda x:
                                                  ','.join([str(y) for y in x]),
                                                  axis=1)
                adata = dutils.filter_cells(adata, 'n_cells', lambda x: x == 9)
            # convert cell lines to mixture id of pure cell line
            else:
                adata.obs['mixture'] = adata.obs['cell_line'].apply(lambda x:
                                    ','.join([str(y) for y in line_to_mixture[x]]))
                adata.obs[lines[0]] = None
                adata.obs[lines[1]] = None
                adata.obs[lines[2]] = None 
                # if you know how to apply for multipel column assignment, lmk
                for x in adata.obs.index.values:
                    cell_type = adata.obs.loc[x, 'cell_line']
                    adata.obs.loc[x, lines] = line_to_mixture[cell_type] 
            adatas[i] = adata
        combined = utils.rbind_adata(adatas)
        sc.pp.filter_genes(combined, min_cells=3)
        # normalize data
        sc.pp.normalize_total(combined)
        # log transform counts to detected highly variable genes
        sc.pp.log1p(combined)
        sc.pp.highly_variable_genes(combined, flavor='seurat')
        # plot highly variable genes
        sc.settings.figdir = snakemake.params['plotdir']
        sc.pl.highly_variable_genes(combined, show=False, save='.svg')
        # subset to highly variable genes
        variable_genes = combined.var.index[
                             np.where(combined.var['highly_variable'])[0]]
        combined = combined[:, variable_genes]
        combined.obs['mixture'] = combined.obs['mixture'].astype('category')
        # write data
        combined.write_csvs(dirname=snakemake.params['outdir'],
                            skip_data=False)
        
                                        
        

