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
        ctrl = snakemake.input['ctrl']
        prtb = snakemake.input['prtb']
        label_col = snakemake.params['label']
        plot_dir = snakemake.params['plotdir']
        out_json = snakemake.output['json']
        
    sc.settings.figdir = plot_dir
    with open(ctrl, 'rb') as f:
        adata = pkl.load(f)
    with open(prtb, 'rb') as f:
        perturbed = pkl.load(f)
    
    performance = main(adata, label_col)
    names = ['controls.svg', 'treated.svg', "combined.svg"]
    adata.obs['Treatment'] = 'Control'
    perturbed.obs['Treatment'] = 'Perturbed'
    combined = utils.rbind_adata([adata, perturbed])
    # add dakota style formatting
    for data, plotfile in zip([adata, perturbed, combined], names):
        fn = os.path.join(plot_dir, 'umap_' + plotfile)
        plot_cells(data, int(performance['n_neighbors']), fn)
        plt.cla()
        plt.clf()
        plt.close()
    with open(out_json, 'w') as f:
        json.dump(performance, f)
