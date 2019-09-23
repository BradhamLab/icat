from scanpy import api as sc
import numpy as np
import pandas as pd

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
        # copy logged.var to combined.var
        adata.var = logged.var
        variable_genes = logged.var.index[
                                np.where(logged.var['highly_variable'])[0]]
        adata = adata[:, variable_genes]
        adata.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)