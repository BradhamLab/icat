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

import colorcet as cc
from icat.src import models, utils
from scanpy import api as sc

try:
    loc = os.path.dirname(os.path.abspath(__file__))
    plt.style.use(os.path.join(loc, 'configs/icat.mplstyle'))
    plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
except:
    pass

def evaluate_icat(control, treated, label_col, icat_kws, plot_dir):
    # adjusted rand no NCFS
    # adjusted rand no ssLouvain
    # adjusted rand with both NCFS and ssLouvain
    # use downstream plotting?
    control.obs['treatment'] = 'Control'
    treated.obs['treatment'] = 'Perturbed'
    combined = sc.AnnData(X=np.vstack([control.X, treated.X]),
                          obs=pd.concat([control.obs, treated.obs],
                                        axis=0,
                                        sort=False).reset_index(drop=True),
                          var=control.var)
    sc.pp.pca(combined)
    try:
        n_neighbors = icat_kws['neighbor_kws']['n_neighbors']
    except KeyError:
        n_neighbors = 15
    sc.tl.pca(combined)
    sc.pp.neighbors(combined, n_neighbors=n_neighbors)
    sc.tl.umap(combined, min_dist=0.0)
    sc.tl.louvain(combined, key_added='louvain')
    icat = models.icat(**icat_kws)
    icat_clustered = icat.cluster(control, treated)
    sc.tl.louvain(icat_clustered, key_added='ncfs-louvain')
    measures = []
    # mutual information, homogeneity score, completeness score
    # rand, fowlkes mallows
    for data, col in zip([combined, icat_clustered, icat_clustered],
                         ['louvain', 'ncfs-louvain', 'sslouvain']):
        print(data)
        split_perf = utils.performance(data, 'Population', col)
        split_perf['method'] = col
        measures.append(split_perf)
        utils.plot_umap(data, col, 'Population')
        plt.savefig(os.path.join(plot_dir, 'umap_' + col + '.svg'))
        plt.cla()
        plt.clf()
        plt.close()
    return measures, icat_clustered


if __name__ == '__main__':
    icat_kws = {'method_kws': {'sigma': 2, 'reg': 3},
                'neighbor_kws': {'n_neighbors': None},
                'cluster_kws': {'resolution': None}}
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        ctrl_X = snakemake.input['ctrl_X']
        ctrl_obs = snakemake.input['ctrl_obs']
        ctrl_var = snakemake.input['ctrl_var']
        prtb_X = snakemake.input['prtb_X']
        prtb_obs = snakemake.input['prtb_obs']
        prtb_var = snakemake.input['prtb_var']
        fit_json = snakemake.input['json']
        name = snakemake.params['name']
        plot_dir = snakemake.params['plotdir']
        out_csv = snakemake.output['csv']
        out_dir = snakemake.params['outdir']
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    control_data = sc.AnnData(X=pd.read_csv(ctrl_X, header=None).values,
                              obs=pd.read_csv(ctrl_obs, index_col=0),
                              var=pd.read_csv(ctrl_var, index_col=0))
    perturb_data = sc.AnnData(X=pd.read_csv(prtb_X, header=None).values,
                              obs=pd.read_csv(prtb_obs, index_col=0),
                              var=pd.read_csv(prtb_var, index_col=0))
    with open(fit_json, 'r') as f:
        fit_data = json.load(f)
    control_data.obs['Population'] = control_data.obs['Population'].astype(str)
    perturb_data.obs['Population'] = control_data.obs['Population'].astype(str)
    icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
    icat_kws['cluster_kws']['resolution'] = fit_data['resolution']
    perf, adata = evaluate_icat(control_data, perturb_data, 'Population',
                                icat_kws, plot_dir)
    adata.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)
    df = pd.DataFrame(perf)
    for key in list(fit_data.keys()):
        if key not in ['n_neighbors', 'resolution']:
            value = fit_data[key]
            df['base.' + key] = value
    parsed_name = utils.parse_sim(name)
    df['Experiment'] = parsed_name['Experiment']
    df['Sim'] = parsed_name['Sim']
    df['Rep'] = parsed_name['Rep']
    df.to_csv(out_csv)
