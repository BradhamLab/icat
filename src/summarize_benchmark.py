import os
import re
import json

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from cycler import cycler
from scanpy import api as sc

from downstream.src.visualization import visualize
from icat.src import plotting, utils

def plot_umaps(adata, plotdir, label, treatment):
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
    visualize.plot_umap(adata, color_col=label)
    plt.savefig(os.path.join(plotdir, 'known_cells_umap.svg'))
    plotting.close_plot()
    # plot umap colored by cluster
    visualize.plot_umap(adata, color_col='Cluster')
    plt.savefig(os.path.join(plotdir, 'cluster_umap.svg'))
    plotting.close_plot()
    # plot umap by 'treatment'
    visualize.plot_umap(adata, color_col=treatment)
    plt.savefig(os.path.join(plotdir, 'treatment_umap.svg'))
    plotting.close_plot()

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        identity = snakemake.params['identity']
        plotdir = snakemake.params['plotdir']
        treatment = snakemake.params['treatment']
        control_id = snakemake.params['controls']
        with open(snakemake.input['fit'], 'r') as f:
            fit = json.load(f)
        n_neighbors = fit['n_neighbors']
        performances = {}
        prefix = os.path.commonpath(snakemake.input['obs']).replace('/', r'\/')
        method_regex = re.compile(r'(?<={}\/)(.*?)(?=\/)'.format(prefix))
        plotted_expectation = False
        for i, each in enumerate(snakemake.input['obs']):
            method = method_regex.search(each).group()
            X = np.loadtxt(snakemake.input['X'][i], delimiter=',')
            obs = pd.read_csv(each, index_col=0)
            adata = sc.AnnData(X=X, obs=obs)
            # standardize cluster column to Cluster
            # rename cluster column to be consistent between methods
            col_name = utils.label_dictionary[method]
            obs.rename(columns={col_name: 'Cluster'}, inplace=True)
            obs['Cluster'] = obs['Cluster'].astype(str)
            # plot umaps
            umap_dir = os.path.join(plotdir, method)
            if not os.path.exists(umap_dir):
                os.makedirs(umap_dir)
            plot_umaps(adata, umap_dir, identity, treatment)
            # plot ternary plots
            performances[method] = utils.performance(obs, identity, 'Cluster')
            # subset data to control and perturbed to asses performance within
            # and out of treatments
            ctrls = obs[obs[treatment] == control_id]
            prtbs = obs[obs[treatment] != control_id]
            for x, subset in zip(['all', 'controls', 'treated'],
                                 [obs, ctrls, prtbs]):
                # ternary plots
                plotting.ternary_plot(subset, 'Cluster', identity)
                plt.savefig(os.path.join(plotdir, '{}_{}_ternary.svg'.\
                                                  format(x, method)))
                plotting.close_plot()
                if not plotted_expectation:
                    plotting.ternary_plot(subset, identity, identity)
                    plt.savefig(os.path.join(plotdir,
                                             '{}_expected_ternary.svg'.\
                                             format(x)))
                    plotting.close_plot()
            plotted_expectation = True
        out = pd.DataFrame(performances)
        out.to_csv(snakemake.output['csv'])
        plt.rcParams['axes.prop_cycle'] =\
            cycler(color=plotting.method_colors)
        plotting.metric_plot(out)
        plt.savefig(snakemake.output['svg'])
        plotting.close_plot()
