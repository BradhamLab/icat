from scanpy import api as sc
import numpy as np
import pandas as pd
from downstream.src.analysis import utils
import os

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        adata = sc.Adata(X=np.loadtxt(snakemake.input['X'], delimiter=','),
                         obs=pd.read_csv(snakemake.input['obs'], index_col=0),
                         var=pd.read_csv(snakemake.input['var'], index_col=0))
        sc.pp.filter_genes(adata, min_cells=3)
        # normalize data
        sc.pp.normalize_total(adata)
        # log transform counts to detected highly variable genes
        logged = adata.copy()
        logged.X = np.log(logged.X + 1)
        sc.pp.highly_variable_genes(logged, flavor='seurat')
        # plot highly variable genes
        sc.settings.figdir = snakemake.params['plotdir']
        sc.pl.highly_variable_genes(logged, show=False, save='.svg')
        # copy logged.var to counts.var
        adata.var = logged.var
        variable_genes = logged.var.index[
                                np.where(logged.var['highly_variable'])[0]]
        adata = adata[:, variable_genes]
        # separate control and perturbed cells
        ctrls = utils.filter_cells(adata, 'stim', lambda x: x == 'ctrl').copy()
        prtbd = utils.filter_cells(adata, 'stim', lambda x: x != 'ctrl').copy()
        # write data to csvs
        ctrls.write_csvs(dirname=os.path.join(snakemake.params['outdir'], 
                                              'controls'),
                        skip_data=False)
        prtbd.write_csvs(dirname=os.path.join(snakemake.params['outdir'], 
                                              'treated'),
                         skip_data=False)