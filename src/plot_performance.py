import itertools
import json
import os

import colorcet as cc
import matplotlib as mpl
import numpy as np
import pandas as pd
from cycler import cycler
from matplotlib import pyplot as plt
from scanpy import api as sc

from downstream.src.visualization import visualize

loc = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(loc, 'configs/figures.mplstyle'))
plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_dark))


method_dictionary = {
    'icat': 'icat',
    'seurat': 'Seurat 2.3 - All',
    'seurat.aligned': 'Seurat 2.3 - Aligned',
    'scanorama': 'scanorama',
    'icat_scan': 'scanorama + icat'
}

label_dictionary = {
    'icat': 'sslouvain',
    'seurat': 'cluster',
    'scanorama': 'scanorama.louvain',
    'icat_scan': 'scanorama.sslouvain'
}

metric_dictionary = {
    'adjusted.mutual.info': 'AMI',
    'adjusted.rand': 'ARI',
    'completeness': 'Completeness',
    'fowlkes.mallows': 'Fowlkes-Mallows',
    'homogeneity': 'Homogeneity'
}

def stacked_barplot(df, label, cluster, xlabel=''):
    """[summary]
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing two columns of sample labels.
    label : str
        Name of column in `df` containing known labels.
    cluster : str
        Name of column in `df` containing predicted labels.
    xlabel : str
        Description of group plotted along the x-axis.
    """
    counts = df.groupby([cluster, label]).size().unstack().fillna(0)
    cluster_totals = counts.sum(axis=0)
    # percentages of cells belonging to each cluster for each known label 
    # (e.g. 87% of cells in known label X are in cluster Y)
    percentages = (counts / cluster_totals).T * 100
    labels = percentages.index.values
    clusters = percentages.columns.values
    xticks = range(percentages.shape[0])
    totals = np.zeros(percentages.shape[0])
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    # plot barplots up to calculated percentage for each cluster in known labels
    for each, color in zip(clusters, colors()):
        new_percentages = percentages.loc[:, each].values
        plt.bar(xticks, new_percentages, color=color['color'],
                width=0.85, bottom=totals, label=each)
        # update cumulate percentage for starting points
        totals += new_percentages
    legend_cols = int(np.ceil(len(clusters) / 20))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
               ncol=legend_cols)
    plt.xticks(xticks, labels, rotation=90)
    plt.xlabel(xlabel, fontsize=24, labelpad=5)
    plt.ylabel("Percentage of Cells", fontsize=24)
    yticks, __ = plt.yticks()
    ylabels = [str(y) + "%" for y in yticks]
    plt.yticks(yticks, ylabels)
    plt.ylim(0, 100)
    plt.tight_layout()

def flip(items, ncol):
    """[summary]

    Taken from here:
    https://stackoverflow.com/questions/10101141/matplotlib-legend-add-items-across-columns-instead-of-down
    
    Parameters
    ----------
    items : [type]
        [description]
    ncol : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def plot_performance(df):
    """
    Plot performance metrics across methods. 
    
    Parameters
    ----------
    df : pd.DataFrame
        Performance data frame where each column is a different performance
        metric, and each row represents a different method.
    """
    indices = np.arange(df.shape[0]).astype(float)
    width = np.min(np.diff(indices)) / df.shape[1]
    indices = indices + indices * width
    colors = cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    starts = indices - width * df.shape[1] / 2
    for method, color in zip(sorted(df.columns), colors()):
        plt.bar(starts + width, df[method], width, color=color['color'],
                label=method_dictionary[method])
        starts = starts + width
    indices = indices + width / 2
    plt.xticks(indices, labels=[metric_dictionary[x] for x in df.index.values])
    plt.title("Method Performance", loc='left')
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    legend_cols = int(np.ceil(df.shape[1] / 2))
    plt.legend(flip(handles, legend_cols), flip(labels, legend_cols),
               loc='upper center', bbox_to_anchor=(0.5, -0.05),
               ncol=legend_cols)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()

def close_plot():
    plt.cla()
    plt.clf() 
    plt.close()

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
        for X, obs, method in zip(snakemake.input['Xs'],
                                  snakemake.input['labels'],
                                  snakemake.params['methods']):
            plotdir = os.path.join(snakemake.params['plotdir'], method)
            if not os.path.exists(plotdir):
                os.makedirs(plotdir)
            X = np.loadtxt(X, delimiter=',')
            obs = pd.read_csv(obs, index_col=0)
            # rename cluster column to be consistent between methods
            col_name = label_dictionary[method]
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
            plt.savefig(os.path.join(plotdir, 'known_cells_umap.svg'))
            close_plot()
            # plot umap colored by cluster
            visualize.plot_umap(adata, color_col='Cluster')
            plt.savefig(os.path.join(plotdir, 'cluster_umap.svg'))
            close_plot()
            # plot umap by 'treatment'
            visualize.plot_umap(adata, color_col=snakemake.params['treatment'])
            plt.savefig(os.path.join(plotdir, 'treatment_umap.svg'))
            close_plot()
            # plot distrubtion of cluster labels between known cell types
            stacked_barplot(adata.obs, snakemake.params['label'], 'Cluster',
                            xlabel=snakemake.params['xlabel'])
            plt.savefig(os.path.join(plotdir, 'cluster_distribution.svg'))
            close_plot()
            # plot_distributin of known cluster labels between clusters
            stacked_barplot(adata.obs, 'Cluster', snakemake.params['label'],
                            xlabel='Cluster')
            plt.savefig(os.path.join(plotdir, 'cell_type_distribution.svg'))
            close_plot()
                  
        # plot performance across metrics
        performance = pd.read_csv(snakemake.input['results'], index_col=0).T
        plot_performance(performance)
        plt.savefig(snakemake.output['metrics'])
        close_plot()
