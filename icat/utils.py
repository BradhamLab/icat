""" 
Utility functions for icat module. 

@author: Dakota Y. Hawkins
@contact: dyh0110@bu.edu 
"""
import collections
import inspect
import sys
import os
import re
import warnings
import logging
import time
from pkg_resources import parse_version

import psutil 
import numpy as np
import pandas as pd
import scanpy as sc 
from scipy import spatial

import igraph as ig

from ncfs import distances 

def ftime(seconds):
    return time.strftime("%Hh:%Mm:%Ss", time.gmtime(seconds))

def set_log():
    # logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    logging.basicConfig(filename='icat.log', level=logging.DEBUG, filemode='w')
    # root = logging.getLogger()
    # root.setLevel(logging.DEBUG)
    # handler = logging.StreamHandler(sys.stdout)
    # handler.setLevel(logging.DEBUG)
    # formatter = logging.Formatter('%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s')
    # handler.setFormatter(formatter)
    # root.addHandler(handler)
    logging.info("System Usage upon startup")
    log_system_usage()

def log_system_usage(msg=None):
    pid = os.getpid()
    py = psutil.Process(pid)
    memory_use = py.memory_info().rss / 2 ** 30 
    if msg is not None:
        logging.info(msg)
    logging.info("Memory usage: {:0.03} GB".format(memory_use))

def close_log(): # this doesn't work
    logging.getLogger().close()

def check_kws(reference_dict, new_dict, name):
    if not isinstance(new_dict, dict):
        raise ValueError("Expected dictionary of keyword arguments for "
                         "`{}`. Received {}.".format(name, type(new_dict)))
    for key, item in new_dict.items():
        if key not in reference_dict.keys():
            raise ValueError("Unsupported keyword argument `{}` for "
                             "{} keywords.".format(key, name))
        reference_dict[key] = item
    return reference_dict


def get_default_kwargs(func, ignore_params=[]):
    params = inspect.signature(func).parameters
    kwargs = {x:params[x].default for x in params if x not in ignore_params}
    return kwargs

def check_string_ids(adata):
    if isinstance(adata.obs.index, pd.RangeIndex):
        warnings.warn("WARNING: Numeric index used for cell ids. " \
                      "Converting to strings.")
        adata.obs.index = adata.obs.index.map(str)
    if isinstance(adata.var.index, pd.RangeIndex):
        warnings.warn("WARNING: Numeric index used for gene ids. " \
                      "Converting to strings.")
        adata.var.index = adata.var.index.map(str)

def check_matching_genes(ref, new):
    """Ensure two AnnData objects have shared genes."""
    return all(ref.var.index == new.var.index)

def subset_cells(adata, vector, value, copy=True):
    """Subset cells by finding indices in `vector` where `value` occurs."""
    idxs = np.where(vector == value)[0]
    if copy:
        return adata[idxs, :].copy()
    return adata[idxs, :]

def check_np_castable(obj, name):
    """Check whether an object is castable to a numpy.ndarray."""
    if not isinstance(obj, np.ndarray):
        try:
            obj = np.array(obj)
        except:
            raise ValueError("Expected numpy.ndarray castable object for "
                            "`{}`. Got {}.".format(name, type(obj)))
    return obj


def get_neighbors(adata, measure='connectivities'):
    # grab connectivities of cells
    if measure not in ['connectivities', 'distances']:
        raise ValueError("`measure` should be either 'connectivities or " \
                         "or 'distances. Received {}.".format(measure))
    if parse_version(sc.__version__) < parse_version("1.5.0"):
        return adata.uns['neighbors'][measure]
    else:
        return adata.obsp[measure]


def distance_matrix(X, metric, kws={}):
    D = np.zeros((X.shape[0], X.shape[0]))
    if isinstance(metric, str):
        try:
            dist_func = distances.supported_distances[metric]
            distances.pdist(X, np.ones(adata.shape[1]),
                            D, dist_func)
        except KeyError:
            D = spatial.distance.pdist(X, metric=metric, **kws)
    else:
        distances.pdist(X, np.ones(adata.shape[1]),
                        D, metric)
    return D



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
    mutables = [True] * len(clusters)
    labels = [None] * len(clusters)
    named_clusters = set([x for x in clusters if not is_none(x)])
    label_to_int = {x: i for i, x in enumerate(sorted(named_clusters))}
    start_label = len(label_to_int)
    for i, x in enumerate(clusters):
        if is_none(x):
            labels[i] = start_label
            start_label += 1
        else:
            labels[i] = label_to_int[x]
            mutables[i] = False
    return (labels, mutables)
