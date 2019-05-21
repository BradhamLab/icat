from scanpy import api as sc
from sklearn import metrics
from icat.src import models
import numpy as np
import pandas as pd
import os
import pickle as pkl
import re

def evaluate_icat(control, treated, label_col, icat_kws):
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
    sc.pp.neighbors(combined, n_neighbors=n_neighbors)
    sc.tl.umap(combined, min_dist=0.0)
    sc.tl.louvain(combined)
    icat = models.icat(**icat_kws)
    icat_clustered = icat.cluster(control, treated)
    sc.tl.louvain(icat_clustered, key_added='ncfs.louvain')
    performance = []
    
    # mutual information, homogeneity score, completeness score
    # rand, fowlkes mallows
    for data, col in zip([combined, icat_clustered, icat_clustered],
                         ['louvain', 'ncfs.louvain', 'sslouvain']):
        known = data.obs[label_col].astype(str)
        pred = data.obs[col].astype(str)
        mi = metrics.adjusted_mutual_info_score(known, pred)
        homog = metrics.homogeneity_score(known, pred)
        comp = metrics.completeness_score(known, pred)
        ar = metrics.adjusted_rand_score(known, pred)
        fm = metrics.fowlkes_mallows_score(known, pred)
        performance.append({'method': col,
                            'adjusted.mutual.info': mi,
                            'homogeneity': homog,
                            'completeness': comp,
                            'adjusted.rand': ar,
                            'fowlkes.mallows': fm})
    return performance


if __name__ == '__main__':
    data_dir = '../data/processed/simulated'
    control_csv = '../data/interim/control_louvain_fits.csv'
    out_csv = '../data/results/icat_performance.csv'
    control_str = 'Controls'
    perturb_str = 'Treated'
    icat_kws = {'method_kws': {'sigma': 2, 'reg': 3},
                'neighbor_kws': {'n_neighbors': None}}
    exp_re = re.compile('^Experiment*[0-9]')
    sim_re = re.compile('Sim*[0-9]')
    rep_re = re.compile('Rep*[0-9]')
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        data_dir = snakemake.params['data']
        control_csv = snakemake.input['csv']
        out_csv = snakemake.output['csv']
        control_str = snakemake.params['ctrl_str']
        perturb_str = snakemake.params['prtb_str']
    fit_data = pd.read_csv(control_csv, index_col=0)
    out_list = []
    for exp in fit_data.index.values[:2]:
        control_fn = os.path.join(data_dir, exp)
        perturb_fn = os.path.join(data_dir, 
                                  exp.replace(control_str, perturb_str))
        with open(control_fn, 'rb') as f:
            control_data = pkl.load(f)
        with open(perturb_fn, 'rb') as f:
            perturb_data = pkl.load(f)
        icat_kws['neighbor_kws']['n_neighbors']\
        = int(fit_data.loc[exp, 'n_neighbors'])
        perf = evaluate_icat(control_data, perturb_data, 'Population', icat_kws)
        df = pd.DataFrame(perf)
        df['Experiment'] = exp_re.search(exp).group()
        df['Sim'] = sim_re.search(exp).group()
        df['Rep'] = sim_re.search(exp).group()
        out_list.append(df)
    out_df = pd.concat(out_list, axis=0).reset_index(drop=True)
    out_df.to_csv(out_csv)
        
    
