import json
import os

import numpy as np
import pandas as pd
from scanpy import api as sc
from matplotlib import pyplot as plt

from downstream.src.visualization import visualize
from icat.src import plotting


label_dictionary = {
    'icat': 'NCFS-SSLouvain',
    'seurat233': 'cluster',
    'scanorama': 'scanorama.louvain',
    'icat_scan': 'scanorama.sslouvain',
    'seurat311': 'seurat_clusters',
    'seurat_icat': 'seurat.sslouvain'
}

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        # load in fit data for n_neighbors
        with open(snakemake.input['fit'], 'r') as f:
            fit = json.load(f)
        n_neighbors = fit['n_neighbors']
        # plot distribution of cluster labels in known 
        plotdir = snakemake.params['plotdir']
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        X = np.loadtxt(snakemake.input['X'], delimiter=',')
        obs = pd.read_csv(snakemake.input['obs'], index_col=0)
        # rename cluster column to be consistent between methods
        col_name = label_dictionary[snakemake.wildcards['method']]
        obs.rename(columns={col_name: 'Cluster'}, inplace=True)
        obs['Cluster'] = obs['Cluster'].astype(str)
        # create AnnData object
        adata = sc.AnnData(X=X, obs=obs)
        # check to see if UMAP projection is already saved in obs data
        # calculate otherwise
        umap_cols = set(['UMAP1', 'UMAP2'])
        if umap_cols.intersection(obs.columns) != umap_cols:
            sc.pp.pca(adata)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
            sc.tl.umap(adata, min_dist=0.0)
            adata.obs['UMAP1'] = adata.obsm['X_umap'][:, 0]
            adata.obs['UMAP2'] = adata.obsm['X_umap'][:, 1]
        # plot umap colored by known cell type
        visualize.plot_umap(adata, color_col=snakemake.params['label'])
        plt.savefig(snakemake.output['known'])
        plotting.close_plot()
        # plot umap colored by cluster
        visualize.plot_umap(adata, color_col='Cluster')
        plt.savefig(snakemake.output['cluster'])
        plotting.close_plot()
        # plot umap by 'treatment'
        visualize.plot_umap(adata, color_col=snakemake.params['treatment'])
        plt.savefig(snakemake.output['treatment'])
        plotting.close_plot()
        # plot distrubtion of cluster labels between known cell types
        plotting.stacked_barplot(adata.obs,
                                 snakemake.params['label'],
                                 'Cluster',
                                 xlabel='Population')
        plt.savefig(snakemake.output['known_bar'])
        plotting.close_plot()
        # plot_distributin of known cell labels between clusters
        plotting.stacked_barplot(adata.obs,
                                 'Cluster',
                                 snakemake.params['label'],
                                 xlabel='Cluster')
        plt.savefig(snakemake.output['cluster_bar'])
        plotting.close_plot()