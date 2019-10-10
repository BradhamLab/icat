import json
import os
import pickle as pkl
import re

import numpy as np
import pandas as pd
import seaborn as sns
from cycler import cycler
from matplotlib import pyplot as plt
from sklearn import metrics

from downstream.src.visualization import visualize
from downstream.src.analysis import utils as dutils

try: 
    import colorcet as cc
except ImportError:
    pass
from icat.src import models, utils
from scanpy import api as sc

try:
    loc = os.path.dirname(os.path.abspath(__file__))
    plt.style.use(os.path.join(loc, 'configs/figures.mplstyle'))
    plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
except:
    pass

def evaluate_icat(control, treated, label_col, icat_kws, plot_dir):
    # adjusted rand no NCFS
    # adjusted rand no ssLouvain
    # adjusted rand with both NCFS and ssLouvain
    combined = utils.rbind_adata([control, treated])
    sc.pp.pca(combined)
    try:
        n_neighbors = icat_kws['neighbor_kws']['n_neighbors']
    except KeyError:
        n_neighbors = 15
    sc.tl.pca(combined)
    sc.pp.neighbors(combined, n_neighbors=n_neighbors)
    sc.tl.umap(combined, min_dist=0.0)
    sc.tl.louvain(combined, key_added='Louvain')
    icat = models.icat(**icat_kws)
    icat_clustered = icat.cluster(control, treated)
    icat_clustered.obs.rename(columns={'sslouvain': 'NCFS-SSLouvain'},
                              inplace=True)
    sc.tl.louvain(icat_clustered, key_added='NCFS-Louvain')
    measures = []
    # mutual information, homogeneity score, completeness score
    # rand, fowlkes mallows
    for data, col in zip([combined, icat_clustered, icat_clustered],
                         ['Louvain', 'NCFS-Louvain', 'NCFS-SSLouvain']):
        split_perf = utils.performance(data, label_col, col)
        split_perf['method'] = col
        measures.append(split_perf)
        visualize.plot_umap(data, color_col=col, shape_col='Population')
        plt.savefig(os.path.join(plot_dir, 'umap_' + col + '.svg'))
        plt.cla()
        plt.clf()
        plt.close()
    return measures, icat_clustered


if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        # load in data
        X = np.loadtxt(snakemake.input['X'], delimiter=',')
        obs = pd.read_csv(snakemake.input['obs'], index_col=0)
        var = pd.read_csv(snakemake.input['var'], index_col=0)
        adata = sc.AnnData(X=X, obs=obs, var=var)
        sc.pp.log1p(adata)
        # subset controls and treated
        ctrls = dutils.filter_cells(adata, snakemake.params['treatment'],
                                    lambda x: x == snakemake.params['control'])\
                                    .copy()
        prtbs = dutils.filter_cells(adata, snakemake.params['treatment'],
                                    lambda x: x != snakemake.params['control'])
        # load in fit data and ncfs ncfs params
        fit_json = snakemake.input['json']
        ncfs_json = snakemake.input['ncfs']
        name = snakemake.params['name']
        plot_dir = snakemake.params['plotdir']
        out_csv = snakemake.output['csv']
        out_dir = snakemake.params['outdir']
        with open(ncfs_json, 'r') as f:
            icat_kws = json.load(f)
        if not os.path.exists(plot_dir):
            os.mkdir(plot_dir)
        with open(fit_json, 'r') as f:
            fit_data = json.load(f)
        # AnnData objects sometimes misbehave with non-string values
        for each in [ctrls, prtbs]:
            each.obs['Population'] = each.obs['Population'].astype(str)
            # each.obs.index = each.obs.index.astype(str)
            # each.var.index = each.var.index.astype(str)
            # each.var.columns = each.var.columns.astype(str)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
        icat_kws['cluster_kws']['resolution'] = fit_data['resolution']
        perf, adata = evaluate_icat(ctrls, prtbs, 'Population',
                                    icat_kws, plot_dir)
        adata.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)
        df = pd.DataFrame(perf)
        for key in list(fit_data.keys()):
            if key not in ['n_neighbors', 'resolution']:
                value = fit_data[key]
                df['base.' + key] = value
        parsed_name = utils.parse_sim(name)
        df['Experiment'] = parsed_name['Experiment']
        df['Perturbation'] = parsed_name['Perturbation']
        df['Sim'] = parsed_name['Sim']
        df['Rep'] = parsed_name['Rep']
        df.to_csv(out_csv)
