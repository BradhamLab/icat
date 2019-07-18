import pandas as pd
from scanpy import api as sc
import json
import numpy as np
import scanorama

def main(adatas, k):
    """
    Integrate scRNAseq datasets using scanorama and cluster with Louvain.

    Parameters
    ----------
    adatas : list, sc.Anndata
        List of single-cell datasets to integrate. Elements should be AnnData
        objects.
    k : int
        Number of neighbors to use during louvain community detection.
    
    Returns
    -------
    [type]
        [description]
    """
    integrated = scanorama.integrate_scanpy(adatas)
    integrated_X = np.vstack(integrated)
    integrated_obs = pd.concat([x.obs for x in adatas], axis=0).\
                        reset_index(drop=True)
    integrated_adata = sc.AnnData(X=integrated_X, obs=integrated_obs)
    sc.pp.pca(integrated_adata)
    sc.pp.neighbors(integrated_adata, n_neighbors=k)
    sc.tl.umap(integrated_adata, min_dist=0.0)
    sc.tl.louvain(integrated_adata, key_added='scanorama.louvain')
    return integrated_adata

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        ctrl = sc.AnnData(X=pd.read_csv(snakemake.input['ctrl_X'],
                                        header=None).values,
                          obs=pd.read_csv(snakemake.input['ctrl_obs'],
                                          index_col=0),
                          var=pd.read_csv(snakemake.input['ctrl_var'],
                                          index_col=0))
        prtb = sc.AnnData(X=pd.read_csv(snakemake.input['prtb_X'],
                                        header=None).values,
                          obs=pd.read_csv(snakemake.input['prtb_obs'],
                                          index_col=0),
                          var=pd.read_csv(snakemake.input['prtb_var'],
                                          index_col=0))
        with open(snakemake.input['json'], 'r') as f:
            fit_data = json.load(f)
        k = fit_data['n_neighbors']
        out = main([ctrl, prtb], k)
        out.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)


