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
        adata = sc.AnnData(X=np.loadtxt(snakemake.input['X'], delimiter=','),
                           obs=pd.read_csv(snakemake.input['obs'], index_col=0),
                           var=pd.read_csv(snakemake.input['var'], index_col=0)) 
        sc.pp.filter_genes(adata, min_cells=3)
        # normalize data
        sc.pp.normalize_total(adata)
        # log transform counts to detected highly variable genes
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor='seurat')
        # plot highly variable genes
        sc.settings.figdir = snakemake.params['plotdir']
        sc.pl.highly_variable_genes(adata, show=False, save='.svg')
        # subset to variable genes
        variable_genes = adata.var.index[
                                np.where(adata.var['highly_variable'])[0]]
        adata = adata[:, variable_genes].copy()
        # write data to csvs
        adata.write_csvs(dirname=snakemake.params['outdir'],
                         skip_data=False)