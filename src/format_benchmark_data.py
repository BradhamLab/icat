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
        line_to_mixture = {lines[0]:'900',
                           lines[1]:'090',
                           lines[2]:'009'}
        bench_regex = re.compile('^(.*?)\.')
        for x, obs in zip(snakemake.input['counts'], snakemake.input['meta']):
            count = pd.read_csv(x, index_col=0)
            obs = pd.read_csv(obs, index_col=0)
            # flag for dataset
            bench = os.path.basename(x[:bench_regex.search(x).end() - 1])
            obs['benchmark'] = os.path.basename(os.path.splitext(x)[0])

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
                adata.obs['cell_line'] = adata.obs[lines].apply(lambda x:
                                                   ''.join([str(y) for y in x]),
                                                   axis=1)
                adata = dutils.filter_cells(adata, 'n_cells', lambda x: x == 9)
            # convert cell lines to mixture id of pure cell line
            else:
                adata.obs['cell_line'] = adata.obs['cell_line'].apply(lambda x:
                                                             line_to_mixture[x])
            adatas[i] = adata
        combined = utils.rbind_adata(adatas)
        # normalize data
        sc.pp.normalize_total(combined)
        for each in adatas:
            bench = each.obs['benchmark'][0]
            out_adata = dutils.filter_cells(combined, 'benchmark',
                                            lambda x: x == bench)
            out_adata.write_csvs(dirname=os.path.join(
                                            snakemake.params['outdir'], bench),
                                 skip_data=False)
        
                                        
        

