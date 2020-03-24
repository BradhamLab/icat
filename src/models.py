"""
Module to identify cell-types across treatments.

@author: Dakota Y. Hawkins
@contact: dyh0110@bu.edu

Example:

.. code-block::python
    from icat import simulate
    data_model = simulate.SingleCellDataset()
    controls = data_model.simulate()
    controls.obs['treatment'] = 'control'
    perturbed = simulate.perturb(controls)
    perturbed.obs['treatment'] = 'perturbed'
    model = icat(ncfs_kws={'reg': 0.5, 'sigma': 3},
                 weight_threshold=1,
                 neighbor_kws={'n_neighbors': 100},
                 sslouvain_kws={'resolution_parameter': 1})
    out = model.cluster(controls, perturbed)
    print(out.obs['sslouvain'])
"""

import inspect
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import sslouvain
from sklearn import preprocessing

from icat.src import utils
from ncfs_expanded import NCFS


class icat():
    """
    Model to identify clusters across treatments.

    Parameters
    ----------
    clustering : str, optional
        Method to use for initial clustering of control cells, by default
        'louvain'. Options are 'louvain' or 'leiden'.
    treatment_col : str, optional
        Column name in oberservation data annotating treatment type between
        datasets, by default 'treatment'.
    ncfs_kws : dict, optional
        Keyword arguments to pass to `ncfs.NCFS` model, by default None and
        default parameters are used. See `ncfs.NCFS` for more information.
    cluster_kws : dict, optional
        Keyword arguments to pass to clustering method, by default None and 
        default parameters are used. See `scanpy.tl.louvain` or
        `scanpy.tl.leiden` for more information.
    cluster_col : str, optional
        Optional column in observation data denoting pre-identified clusters.
        By default None, and clustering will be performed following
        `clustering` and `cluster_kws`. 
    weight_threshold : float, optional
        Heuristic to determine number of informative genes returned by NCFS.
        Does not affect performance, but does help inform reasonable NCFS 
        hyperparamters. Default is 1. 
    neighbor_kws : dict, optional
        Keyword arguments for identifying neighbors. By default None, and 
        default parameters are used. See `scanpy.pp.neighbors` for more
        information.
    sslouvain_kws : dict, optional
        Keyword arguments for `sslouvain`. By default None, and default
        parameters will be used. See `sslouvain.find_partition` for more
        information.
    pca_kws : dict, optional
        Keyword arguments for performing principle component analysis. B default
        None, and default parameters will be used. See `sslouvain.tl.pca` for
        more information.
    """

    def __init__(self, clustering='louvain', treatment_col='treatment',
                 ncfs_kws=None, cluster_kws=None, cluster_col=None,
                 weight_threshold=1.0, neighbor_kws=None,
                 sslouvain_kws=None, pca_kws=None):
        self.clustering = clustering
        self.treatment_col = treatment_col
        self.ncfs_kws = ncfs_kws
        self.cluster_kws = cluster_kws
        self.cluster_col = cluster_col
        self.weight_threshold = weight_threshold
        self.neighbor_kws = neighbor_kws
        self.sslouvain_kws = sslouvain_kws
        self.pca_kws = pca_kws
    
    @property
    def clustering(self):
        return self._clustering

    @clustering.setter
    def clustering(self, value):
        if not isinstance(value, str):
            raise ValueError("Expected string for `clustering` parameter. "
                             "Received {}".format(value))
        value = value.lower()
        if value not in ['louvain', 'leiden']:
            raise ValueError("Unsupported method: {}".format(value))
        self._clustering = value

    @property
    def treatment_col(self):
        """
        Column in control and perturbed datasets annotated treatment.
        """
        return self._treatment_col

    @treatment_col.setter
    def treatment_col(self, value):
        self._treatment_col = value

    @property
    def ncfs_kws(self):
        """
        Keyword arguments to pass to NCFS. See NCFS documentation for me information.
        """
        return self._ncfs_kws
    
    @ncfs_kws.setter
    def ncfs_kws(self, value):
        default_kws = utils.get_default_kwargs(NCFS.NCFS, ['self'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'ncfs')
        else:
            value = default_kws
        self._ncfs_kws = value
    
    @property
    def weight_threshold(self):
        """ Weight threshold for genes to be considered informative."""
        return self._weight_threshold
    
    @weight_threshold.setter
    def weight_threshold(self, value):
        elif not isinstance(value, (float, int, np.float, np.integer)):
            raise ValueError("Expected numerical value for `weight_threshold`."
                             "Received: {}".format(value))
        self._weight_threshold = value

    @property
    def neighbor_kws(self):
        """Keyword arguments for neighbor determination."""
        return self._neighbor_kws

    @neighbor_kws.setter
    def neighbor_kws(self, value):
        default_kws = utils.get_default_kwargs(sc.pp.neighbors,
                                               ['adata', 'copy'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'neighbor_kws')
        else:
            value = default_kws
        self._neighbor_kws = value
        
    @property
    def cluster_kws(self):
        """Keyword arguments to pass to clustering algorithm."""
        return self._cluster_kws
    
    @cluster_kws.setter
    def cluster_kws(self, value):
        default_kws = {'louvain': utils.get_default_kwargs(sc.tl.louvain,
                                                     ['adata', 'copy'])}
        if 'leiden' in dir(sc.tl):
            default_kws['leiden'] = utils.get_default_kwargs(sc.tl.leiden,
                                                             ['adata', 'copy'])
        else:
            if self.clustering == 'leiden':
                raise ValueError('The version of Scanpy you are running does '\
                                 'not support Leiden clustering. If you wish '\
                                 'to use Leiden clustering, try upgrading your'\
                                 'version of Scanpy.')

        if value is not None:
            value = utils.check_kws(default_kws[self.clustering],
                                    value, 'cluster.' + self.clustering)
        else:
            value = default_kws[self.clustering]
        self._cluster_kws = value

    @property
    def cluster_col(self):
        """Optional column identifying pre-computed clusters."""
        return self._cluster_col
    
    @cluster_col.setter
    def cluster_col(self, value):
        if not isinstance(value, (str, int)) and value is not None:
            raise ValueError('Expected integer index or string name for '+\
                             '`cluster_col`')
        self._cluster_col = value

    @property
    def sslouvain_kws(self):
        """Keyword arguments for sslouvain."""
        return self._sslouvain_kws
    
    @sslouvain_kws.setter
    def sslouvain_kws(self, value):
        default_kws = utils.get_default_kwargs(sslouvain.find_partition, '')
        if value is not None:
            value = utils.check_kws(default_kws, value, 'sslouvain_kws')
        else:
            value = default_kws
        self._sslouvain_kws = value

    @property
    def pca_kws(self):
        """Keyword arguments for PCA."""
        return self._pca_kws

    @pca_kws.setter
    def pca_kws(self, value):
        default_kws = utils.get_default_kwargs(sc.pp.pca, ['data', 'copy'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'pca_kws')
        else:
            value = default_kws
        self._pca_kws = value

    def cluster(self, controls, perturbed):
        """
        Cluster cells in control and experimental conditions.
        
        Parameters
        ----------
        controls : sc.AnnData
            Annotated dataframe of control cells.
        perturbed : sc.AnnData
            Annotated dataframe of treated/perturbed cells.
        
        Returns
        -------
        sc.AnnData
            Annotated dataframe of combined cells in NCFS space. 
        """
        if not isinstance(controls, sc.AnnData):
            raise ValueError("Expected AnnData object for `controls`.")
        if self.treatment_col not in controls.obs.columns:
            controls.obs[self.treatment_col] = 'Control'
        if not isinstance(perturbed, sc.AnnData):
            if isinstance(perturbed, list):
                if not all([isinstance(x, sc.AnnData) for x in perturbed]):
                    raise ValueError("Expected all perturbed datasets to be "
                                     "sc.AnnData objects.")
                if not all([utils.check_matching_genes(controls, x)\
                for x in perturbed]):
                    raise ValueError("Gene columns do not match between control"
                                     " and perturbed cells.")
                perturbed = perturbed[0].concatenate(*perturbed[1:])
            else:
                raise ValueError("Unexpected input type for `perturbed`: "
                                 "{}. Expected list of sc.AnnData objects or "
                                 "a single sc.AnnData object".\
                                 format(type(perturbed)))
        if utils.check_matching_genes(controls, perturbed):
            raise ValueError("Gene columns do not match between control and"
                                " perturbed cells.")
        if self.treatment_col not in perturbed.obs.columns:
            raise ValueError("Expected {} column in perturbed data.".format(
                                self.treatment_col))
        
        # change numeric indices to strings
        for each in [controls, perturbed]:
            if isinstance(each.obs.index, pd.RangeIndex):
                print("WARNING: Numeric index used for cell ids. "
                      "Converting to strings.")
                each.obs.index = each.obs.index.map(str)
            if isinstance(each.var.index, pd.RangeIndex):
                print("WARNING: Numeric index used for gene ids. "
                      "Converting to strings.")
                each.var.index = each.var.index.map(str)
        # no previous clustering provided, cluster using louvain or leiden
        if self.cluster_col is None:
            if self.neighbor_kws['use_rep'] != 'X':
                sc.pp.pca(controls, **self.pca_kws)
            sc.pp.neighbors(controls, **self.neighbor_kws)
            sc.tl.umap(controls, min_dist=0.0)
            if self.clustering == 'louvain':
                sc.tl.louvain(controls, **self.cluster_kws)
                self.cluster_col = 'louvain'
            elif self.clustering == 'leiden':
                sc.tl.leiden(controls, **self.cluster_kws)
                self.cluster_col = 'leiden'
        else:
            if self.cluster_col not in controls.obs.columns:
                raise ValueError(f"`cluster_col` - {self.cluster_col} not found"
                                  " in control data.")
            if np.any([pd.isnull(x) for x in controls.obs[self.cluster_col]]):
                raise ValueError("Expected labelled cells by passing "\
                                 "`cluster_col`={}. ".format(self.cluster_col) +
                                 "Received atleast one unannotated cell.")
        # combine control and perturbed data
        combined = controls.concatenate(perturbed, join='outer')
        # fit scale function to scale features between 0 and 1
        scaler = preprocessing.MinMaxScaler()
        scaler.fit(controls.X)
        # instantiate ncfs model
        model = NCFS.NCFS(**self.ncfs_kws)
        fit_X = scaler.transform(controls.X).astype(np.float64)
        
        # fit gene weights using control dataset
        model.fit(fit_X, np.array(controls.obs[self.cluster_col].values),
                  sample_weights='balanced')
        # scale combined matrix between 0 and 1 and apply learned weights
        # across gene expression matrix
        combined.X = model.transform(scaler.transform(combined.X))
        # save genes weights
        combined.var['ncfs.weights'] = model.coef_
        combined.var['informative'] = combined.var['ncfs.weights'].apply(lambda x: x > self.weight_threshold)
        n_clusters = len(controls.obs[self.cluster_col].unique())
        informative = combined.var['informative'].sum()
        if sum(model.coef_ > self.weight_threshold) < n_clusters:
            warnings.warn("Number of informative genes less "
                          "than the number of identified control clusters: "
                          f"informative genes: {informative}, "
                          f"number of clusters:  {n_clusters}. Consider "
                          "increasing `sigma` or decreasing `reg` for better "
                          "performance.")
        # create neighbor graph for control+perturbed combined data
        self.neighbor_kws['use_rep'] = 'X'
        sc.pp.neighbors(combined, **self.neighbor_kws)
        sc.tl.umap(combined, min_dist=0.0)
        # grab connectivities of cells
        g = utils.igraph_from_adjacency(combined.uns['neighbors']['connectivities'])
        # instantiate semi-supervised Louvain model
        try:
            resolution = self.sslouvain_kws['resolution_parameter']
        except KeyError:
            resolution = 1.0
        y_, mutables = utils.format_labels(combined.obs[self.cluster_col])
        part = sslouvain.find_partition(g,
                                        sslouvain.RBConfigurationVertexPartition,
                                        initial_membership=y_,
                                        mutable_nodes=mutables,
                                        resolution_parameter=resolution)
        # store new cluster labels in cell metadata
        combined.obs['sslouvain'] = part.membership
        return combined

