import scanpy.api as sc
import matplotlib.pyplot as plt
from cycler import cycler
import pickle as pkl
from icat.src import utils
import json
import pandas as pd
import numpy as np
import os
import colorcet as cc
from downstream.src.visualization import visualize
from downstream.src.analysis import utils as dutils

try:
    loc = os.path.dirname(os.path.abspath(__file__))
    plt.style.use(os.path.join(loc, 'configs/figures.mplstyle'))
    plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
except:
    pass

def main(adata, label_col):
    sc.pp.pca(adata)
    performance = dict()
    n_max = int(0.60 * adata.shape[0])
    n_min = 3
    score = -np.inf
    resolutions = [1]
    for n in range(n_min, n_max + 1):
        for r in resolutions:
            sc.pp.neighbors(adata, n_neighbors=n)
            sc.tl.umap(adata, min_dist=0.0)
            sc.tl.louvain(adata, resolution=r, key_added='louvain')
            measures = utils.performance(adata, label_col, 'louvain')
            if score < measures['adjusted.rand']:
                performance = measures.copy()
                performance['n_neighbors'] = n
                performance['resolution'] = r
                score = measures['adjusted.rand']
    return performance


def plot_cells(adata, n, fn, color, shape):
    adata.obs[color] = adata.obs[color].astype(str).astype('category')
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=n)
    sc.tl.umap(adata, min_dist=0)
    adata.obs['UMAP1'] = adata.obsm['X_umap'][:, 0]
    adata.obs['UMAP2'] = adata.obsm['X_umap'][:, 1]
    visualize.plot_umap(adata, color_col=color, shape_col=shape)
    plt.savefig(fn)
    plt.cla()
    plt.clf()
    plt.close()

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        X = np.loadtxt(snakemake.input['X'], delimiter=',')
        print(X.shape)
        obs = pd.read_csv(snakemake.input['obs'], index_col=0)
        print(obs.shape)
        var = pd.read_csv(snakemake.input['var'], index_col=0)
        print(var.shape)
        adata = sc.AnnData(X=X, obs=obs, var=var)
        separated = {y: dutils.filter_cells(adata,
                                            snakemake.params['treatment'],
                                            lambda x: x==y).copy()\
                    for y in adata.obs[snakemake.params['treatment']].unique()}
        separated['combined'] = adata
        performance = main(separated[snakemake.params['control']],
                           snakemake.params['label'])
        fit_data = main(separated[snakemake.params['control']],
                        snakemake.params['label'])
        if not os.path.exists(snakemake.params['plotdir']):
            os.makedirs(snakemake.params['plotdir'])
        # add dakota style formatting
        for treatment, data in separated.items():
            plotfile = os.path.join(snakemake.params['plotdir'],
                                    "umap_{}.svg".format(treatment))
            plot_cells(data, int(performance['n_neighbors']), plotfile,
                       snakemake.params['label'], snakemake.params['treatment'])
        with open(snakemake.output['json'], 'w') as f:
            json.dump(fit_data, f)
