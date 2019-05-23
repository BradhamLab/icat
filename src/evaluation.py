from scanpy import api as sc
from sklearn import metrics
from icat.src import models
import numpy as np
import pandas as pd
import os
import pickle as pkl
import re
import json
from icat.src import utils
from matplotlib import pyplot as plt
import seaborn as sns
from cycler import cycler
import colorcet as cc

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
    return measures


if __name__ == '__main__':
    icat_kws = {'method_kws': {'sigma': 2, 'reg': 3},
                'neighbor_kws': {'n_neighbors': None},
                'cluster_kws': {'resolution': None}}
    exp_re = re.compile('^Experiment*[0-9]')
    sim_re = re.compile('Sim*[0-9]')
    rep_re = re.compile('Rep*[0-9]')
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        ctrl = snakemake.input['ctrl']
        prtb = snakemake.input['prtb']
        fit_json = snakemake.input['json']
        name = snakemake.params['name']
        plot_dir = snakemake.params['plotdir']
        out_csv = snakemake.output['csv']
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    with open(ctrl, 'rb') as f:
        control_data = pkl.load(f)
    with open(prtb, 'rb') as f:
        perturb_data = pkl.load(f)
    with open(fit_json, 'r') as f:
        fit_data = json.load(f)
    control_data.obs['Population'] = control_data.obs['Population'].astype(str)
    perturb_data.obs['Population'] = control_data.obs['Population'].astype(str)
    icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
    icat_kws['cluster_kws']['resolution'] = fit_data['resolution']
    perf = evaluate_icat(control_data, perturb_data, 'Population', icat_kws,
                         plot_dir)
    df = pd.DataFrame(perf)
    for key in list(fit_data.keys()):
        if key not in ['n_neighbors', 'resolution']:
            value = fit_data[key]
            df['base.' + key] = value
    df['Experiment'] = exp_re.search(name).group()
    df['Sim'] = sim_re.search(name).group()
    df['Rep'] = rep_re.search(name).group()
    df.to_csv(out_csv)
        
    
