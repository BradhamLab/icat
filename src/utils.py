import collections
import inspect
import os
import re
import warnings

import numpy as np
import pandas as pd
import scanpy.api as sc
import seaborn as sns
from cycler import cycler
from matplotlib import pyplot as plt
from sklearn import metrics
import igraph as ig


# try:
#     import colorcet as cc
#     loc = os.path.dirname(os.path.abspath(__file__))
#     plt.style.use(os.path.join(loc, 'configs/icat.mplstyle'))
#     plt.rc('axes', prop_cycle=cycler('color', cc.glasbey_light))
# except ImportError:
#     pass

label_dictionary = {
    'icat': 'sslouvain',
    'seurat233': 'cluster',
    'scanorama': 'scanorama.louvain',
    'icat_scan': 'scanorama.sslouvain',
    'seurat311': 'seurat_clusters',
    'seurat_icat': 'seurat.sslouvain'
}


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
                                        sort=False),
                          var=pd.concat([each.var[col] for each, col\
                                        in zip(adata_lists, df_cols)], axis=1))
    if not combined.obs.index.is_unique:
        print("WARNING: Datasets contain overlapping" 
              " cell ids. Index will be reset.")
        combined.obs.index = ['cell-{}'.format(i + 1)\
                              for i in range(combined.shape[0])]
    return combined


def performance(adata, true_col, pred_col):
    """Measure the performance of a clustering partition."""
    if isinstance(adata, sc.AnnData):
        obs = adata.obs
    elif isinstance(adata, pd.DataFrame):
        obs = adata
    else:
        raise TypeError("Unsupported type: {}".format(type(adata)))
    known = obs[true_col].values.astype(str)
    pred = obs[pred_col].values.astype(str)
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

def igraph_from_adjacency(adjacency, directed=None):
    """
    Get igraph graph from adjacency matrix.
    
    Parameters
    ----------

    adjacency : numpy.ndarray, scipy.sparse.csr
        Adjacency matrix where non-zero entries represent connections between 
        samples.
    directed: boolean, optional

    Returns
    -------
    igraph.graph
        Nearest neighbor graph generated from adjacency matrix.

    References
    ----------
    Taken from: https://github.com/theislab/scanpy/blob/28498953092dc7cbecd0bd67380b1b060367d639/scanpy/_utils.py#L170
    """
    import igraph as ig
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        warnings.warn(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g

def is_none(x):
    """Check whether a value is a null value."""
    if isinstance(x, float):
        return np.isnan(x)
    return x is None

def format_labels(clusters):
    """
    Format cluster labels for sslouvain.

    Parameters
    ----------
    clusters : iterable
        List of cluster labels where None or np.nan represent previously
        unlabeled samples.
    
    Returns
    -------
    (list, list):
        labels : list
            List of labels where previously unlabeled samples are given 
            their own labels
        mutables : list
            List of samples indicating which samples were previously labelled.
    """
    if isinstance(clusters, pd.Series):
        clusters = clusters.values
    elif isinstance(clusters, list):
        clusters = np.array(clusters)
    elif not isinstance(clusters, np.ndarray):
        warnings.warn(f"Unsupport type {type(clusters)} for `clusters`")
    clusters = clusters.astype(float)
    mutables = [True] * len(clusters)
    labels = [None] * len(clusters)
    start_label = int(np.nanmax(clusters) + 1)
    for i, x in enumerate(clusters):
        if is_none(x):
            labels[i] = start_label
            start_label += 1
        else:
            labels[i] = int(x)
            mutables[i] = False
    return (labels, mutables)

def plot_umap(adata, color, shape, ax=None):
    if ax is None:
        __, ax = plt.subplots(figsize=(10, 8))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
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
    exp_re = re.compile('^Experiment.[0-9a-z]|Experiment*[0-9]')
    sim_re = re.compile('Sim*[0-9]')
    rep_re = re.compile('Rep*[0-9]')
    pert_re = re.compile('Perturbation*[0-9]')
    out = {}
    # scan name for experiment, perturbation, sim, and rep ids
    # if no match, return empty string
    for x, regex in zip(['Experiment', 'Perturbation', 'Sim', 'Rep'],
                        [exp_re, pert_re, sim_re, rep_re]):
        match = regex.search(name)
        if match is not None:
            out[x] = match.group()
        else:
            out[x] = ''
        if x == 'Experiment' and out[x] == '':
            with open('parse.log', 'a') as f:
                f.write('name: {}\n'.format(name))
    return out
