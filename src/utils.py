import collections
import inspect
import os
import re

import numpy as np
import pandas as pd
import scanpy.api as sc
import seaborn as sns
from cycler import cycler
from matplotlib import pyplot as plt
from sklearn import metrics

import colorcet as cc

try:
    loc = os.path.dirname(os.path.abspath(__file__))
    plt.style.use(os.path.join(loc, 'configs/icat.mplstyle'))
    plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
except:
    pass

def check_kws(reference_dict, new_dict, name):
    if not isinstance(new_dict, dict):
        raise ValueError("Expected dictionary of keyword arguments for "
                         "`{}`. Received {}.".format(name, type(new_dict)))
    for key, item in new_dict.items():
        if key not in reference_dict.keys():
            raise ValueError("Unsupported keyword argument `{}` for "
                             "{} keywords.".format(key, name))
        new_dict[key] = item
    return new_dict


def check_matching_genes(ref, new):
    return set(ref.var.index.values).difference(new.var.index.values) == 0


def get_default_kwargs(func, ignore_params=[]):
    params = inspect.signature(func).parameters
    kwargs = {x:params[x].default for x in params if x not in ignore_params}
    return kwargs


def check_np_castable(obj, name):
    """Check whether an object is castable to a numpy.ndarray."""
    if not isinstance(obj, np.ndarray):
        try:
            obj = np.array(obj)
        except:
            raise ValueError("Expected numpy.ndarray castable object for "
                            "`{}`. Got {}.".format(obj, type(obj)))
    return obj


def __evaluate_key(key, sep):
    if not isinstance(key, str):
        raise ValueError("Keys must be strings for easy concatentation.")
    if sep in key:
        raise ValueError("Cannot have `{}` in dictionary keys.".format(sep))


def flatten_dict(d, parent_key='', sep='.'):
    """
    Flatten a dictionary containing nested dictionaries.
    
    Parameters
    ----------
    d : dict
        Dictionary to flatten
    parent_key : str, optional
        Key in parent dictionary pointing to `d`. The default is '', which
        assumes `d` is the highest level nested dictionary.
    sep : str, optional
        String value to separate child and parent keys. The default is '.',
        which will place a '.' between each key. All parent and child keys
        will be assessed to ensure they do not contain a `sep` character;
        therefore, `sep` should be set to a delimiter not present in current
        keys.
    
    Returns
    -------
    dict
        Flattened dictionary with parent and child keys separted by `sep`.

    References
    ----------

    Taken shamelessly from here:
        https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
    """

    items = []
    for k, v in d.items():
        __evaluate_key(k, sep)
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def rbind_adata(adata_lists):
    df_cols = [adata_lists[0].var.columns] + [None]*(len(adata_lists) - 1)
    for i, adata in enumerate(adata_lists[1:]):
        df_cols[i + 1] = adata.var.columns.difference(df_cols[i])
    combined = sc.AnnData(X=np.vstack([each.X for each in adata_lists]),
                          obs=pd.concat([each.obs for each in adata_lists],
                                        axis=0,
                                        sort=False).reset_index(drop=True),
                          var=pd.concat([each.var[col] for each, col\
                                        in zip(adata_lists, df_cols)], axis=1))
    return combined


def performance(adata, true_col, pred_col):
    """Measure the performance of a clustering partition."""
    known = adata.obs[true_col].values.astype(str)
    pred = adata.obs[pred_col].values.astype(str)
    mi = metrics.adjusted_mutual_info_score(known, pred, 
                                            average_method='arithmetic')
    homog = metrics.homogeneity_score(known, pred)
    comp = metrics.completeness_score(known, pred)
    ar = metrics.adjusted_rand_score(known, pred)
    fm = metrics.fowlkes_mallows_score(known, pred)
    measures = {'adjusted.mutual.info': mi,
                'homogeneity': homog,
                'completeness': comp,
                'adjusted.rand': ar,
                'fowlkes.mallows': fm}
    return measures


def plot_umap(adata, color, shape, ax=None):
    if ax is None:
        __, ax = plt.subplots(figsize=(10, 8))
    colors =plt.rcParams['axes.prop_cycle'].by_key()['color']
    pallete = {}
    adata.obs[color] = adata.obs[color].astype(str)
    for i, each in enumerate(adata.obs[color].unique()):
        pallete[each] = colors[i]
    figure = sns.scatterplot(x=adata.obsm['X_umap'][:, 0],
                             y=adata.obsm['X_umap'][:, 1],
                             hue=adata.obs[color],
                             palette=pallete,
                             style=adata.obs[shape].astype(str),
                             ax=ax,
                             s=150)
    legend = ax.get_legend()
    for i, handle in enumerate(legend.legendHandles):
        if handle.get_label() == shape:
            legend.legendHandles[i].set_facecolor('white')
            legend.legendHandles[i].set_color('white')
            legend.legendHandles[i].set_edgecolor('white')

    return figure

def parse_sim(name):
    exp_re = re.compile('^Experiment*[0-9]')
    sim_re = re.compile('Sim*[0-9]')
    rep_re = re.compile('Rep*[0-9]')
    
    out = {'Experiment': exp_re.search(name).group(),
           'Sim': sim_re.search(name).group(),
           'Rep': rep_re.search(name).group()}
    return out
