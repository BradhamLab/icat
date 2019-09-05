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

try:
    loc = os.path.dirname(os.path.abspath(__file__))
    plt.style.use(os.path.join(loc, 'configs/icat.mplstyle'))
    plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
except:
    pass

def main(andata, label_col):
    sc.pp.pca(andata)
    performance = dict()
    n_max = int(0.60 * andata.shape[0])
    n_min = 3
    score = -np.inf
    resolutions = [1]
    for n in range(n_min, n_max + 1):
        for r in resolutions:
            sc.pp.neighbors(andata, n_neighbors=n)
            sc.tl.umap(andata, min_dist=0.0)
            sc.tl.louvain(andata, resolution=r, key_added='louvain')
            measures = utils.performance(andata, label_col, 'louvain')
            if score < measures['adjusted.rand']:
                performance = measures.copy()
                performance['n_neighbors'] = n
                performance['resolution'] = r
    return performance


def plot_cells(andata, n, fn):
    andata.obs['Population'] = andata.obs['Population'].astype('category')
    sc.pp.pca(andata)
    sc.pp.neighbors(andata, n_neighbors=n)
    sc.tl.umap(andata, min_dist=0)
    utils.plot_umap(andata, color='Population', shape='Treatment')
    plt.savefig(fn)
    plt.cla()
    plt.clf()
    plt.close()

if __name__ == '__main__':
    # ctrl = '../data/processed/simulated/Experiment1Sim1Rep1-Controls.pkl'
    # prtb = '../data/processed/simulated/Experiment1Sim1Rep1-Treated.pkl'
    # label_col = 'Population'
    # out = '../data/interim/Experiment1Sim1Rep1-Controls_fits.json'
    # out = '../figures/simulated/Experiment1Sim1Rep1/'
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        ctrl_X = snakemake.input['ctrl_X']
        ctrl_obs = snakemake.input['ctrl_obs']
        try:
            ctrl_var = snakemake.input['ctrl_var']
        except AttributeError:
            ctrl_var = None
        prtb_X = snakemake.input['prtb_X']
        prtb_obs = snakemake.input['prtb_obs']
        try:
            prtb_var = snakemake.input['prtb_var']
        except AttributeError:
            prtb_var = None
        label_col = snakemake.params['label']
        plot_dir = snakemake.params['plotdir']
        out_json = snakemake.output['json']
        
    sc.settings.figdir = plot_dir
    if ctrl_var is not None:
        ctrl = sc.AnnData(X=pd.read_csv(ctrl_X, header=None).values,
                          obs=pd.read_csv(ctrl_obs, index_col=0),
                          var=pd.read_csv(ctrl_var, index_col=0))
    else:
        ctrl = sc.AnnData(X=pd.read_csv(ctrl_X, header=None).values,
                          obs=pd.read_csv(ctrl_obs, index_col=0))
    if isinstance(prtb_X, list):
        adatas = []
        for x, obs in zip(sorted(prtb_X), sorted(prtb_obs)):
            adata = sc.AnnData(X=np.loadtxt(x, delimiter=','),
                               obs=pd.read_csv(obs, index_col=0))
            adatas.append(adata)
        prtb = utils.rbind_adata(adatas)
        ctrl.obs['Mixture'] = 'No'
        prtb.obs['Mixture'] = 'Yes'
    
    else:
        prtb = sc.AnnData(X=pd.read_csv(prtb_X, header=None).values,
                          obs=pd.read_csv(prtb_obs, index_col=0),
                          var=pd.read_csv(prtb_var, index_col=0))
        ctrl.obs['Treatment'] = 'Control'
        prtb.obs['Treatment'] = 'Perturbed'
    
    performance = main(ctrl, label_col)
    names = [snakemake.output['ctrl_svg'], snakemake.output['prtb_svg'],
             snakemake.output['comb_svg']]
    combined = utils.rbind_adata([ctrl, prtb])
    # add dakota style formatting
    for data, plotfile in zip([ctrl, prtb, combined], names):
        fn = os.path.join(plot_dir, plotfile)
        plot_cells(data, int(performance['n_neighbors']), fn)
        plt.cla()
        plt.clf()
        plt.close()
    with open(out_json, 'w') as f:
        json.dump(performance, f)
