import pandas as pd
from scanpy import api as sc
import json
import numpy as np
import scanorama
from icat.src import utils


def run_scanorama(adatas, k):
    integrated_X = scanorama.integrate_scanpy(adatas)
    integrated = []
    for i, each in enumerate(integrated_X):
        adata = sc.AnnData(X=each, obs=adatas[i].obs)
        integrated.append(adata)
    data = utils.rbind_adata(integrated)
    sc.pp.pca(data)
    sc.pp.neighbors(data, n_neighbors=k)
    sc.tl.umap(data, min_dist=0.0)
    sc.tl.louvain(data, key_added='scanorama.louvain')
    return data

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        all_X = sorted(snakemake.input['X'])
        all_obs = sorted(snakemake.input['obs'])
        all_var = sorted(snakemake.input['var'])
        adatas = []
        for i, each in enumerate(all_X):
            adata = sc.AnnData(X=np.loadtxt(each, delimiter=','),
                               obs=pd.read_csv(all_obs[i], index_col=0),
                               var=pd.read_csv(all_var[i], index_col=0))
            if snakemake.params['control_id'] in each:
                adatas = [adata] + adatas
            else:
                adatas.append(adata)
        with open(snakemake.input['json'], 'r') as f:
            fit_data = json.load(f)
        k = fit_data['n_neighbors']
        out = run_scanorama(adatas, k)
        out.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)