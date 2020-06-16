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
import logging
import time
from pkg_resources import parse_version

import numpy as np
import pandas as pd
import scanpy as sc
import sslouvain
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

import ncfs

from . import utils


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
    reference : str, optional
        Which data partitions to use for when learning gene weights. Default is 
        'all' and all partions by treatment type will be used (as defined by 
        values in `treatment_col`). Acceptable values are 'all' and 'controls'.
    sub_sample : bool, optional
        Whether to perform NCFS training on a sample of the data. Default is
        False, and training will be performed on the entire dataset.
    train_size : float, optional.
        Proportion of data to train on if `sub_sample==True`. Values should be
        between 0 and 1. Default is 0.1, and training will be performed on 10%
        of the data.
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
                 reference='all', sub_sample=False, train_size=0.1,
                 ncfs_kws=None, cluster_kws=None, cluster_col=None,
                 weight_threshold=1.0, neighbor_kws=None,
                 sslouvain_kws=None, pca_kws=None):
        self.clustering = clustering
        self.treatment_col = treatment_col
        self.reference = reference
        self.sub_sample = sub_sample
        self.train_size = train_size
        self.ncfs_kws = ncfs_kws
        self.cluster_kws = cluster_kws
        self.cluster_col = cluster_col
        self.weight_threshold = weight_threshold
        self.neighbor_kws = neighbor_kws
        self.sslouvain_kws = sslouvain_kws
        self.pca_kws = pca_kws
        # utils.set_log()
    
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
    def reference(self):
        """Reference data set. Acceptable values 'all' or 'controls'."""
        return self._reference

    @reference.setter
    def reference(self, value):
        if value not in ['all', 'controls']:
            raise ValueError("`reference` value should be either 'all' or 'controls'")
        self._reference = value

    @property
    def sub_sample(self):
        "Whether to perform NCFS training on a sample of the data."
        return self._sub_sample
    
    @sub_sample.setter
    def sub_sample(self, value):
        if not isinstance(value, bool):
            raise ValueError('Expected bool for `sub_sample`')
        self._sub_sample = value

    @property
    def train_size(self):
        "Proportion of data to train on. Default is 0.1"
        return self._train_size
    
    @train_size.setter
    def train_size(self, value):
        if not 0 < value < 1:
            raise ValueError("Expected value between 0 and 1 for `train_size`.")
        self._train_size = value

    @property
    def ncfs_kws(self):
        """
        Keyword arguments to pass to NCFS. See NCFS documentation for me information.
        """
        return self._ncfs_kws
    
    @ncfs_kws.setter
    def ncfs_kws(self, value):
        default_kws = utils.get_default_kwargs(ncfs.NCFS, ['self'])
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
        if not isinstance(value, (float, int, np.float, np.integer)):
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
        default_kws['metric'] = ncfs.distances.phi_s
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

    def cluster(self, controls, perturbed, verbose=False):
        """
        Cluster cells in control and experimental conditions.
        
        Parameters
        ----------
        controls : sc.AnnData
            Annotated dataframe of control cells.
        perturbed : sc.AnnData
            Annotated dataframe of treated/perturbed cells.
        verbose : boolean, optional
            Whether to print progress through clustering step. 
        
        Returns
        -------
        sc.AnnData
            Annotated dataframe of combined cells in NCFS space. 
        """
        if not isinstance(controls, sc.AnnData):
            raise ValueError("Expected AnnData object for `controls`.")
        if self.treatment_col not in controls.obs.columns:
            controls.obs[self.treatment_col] = 'Control'
        # change numeric indices to strings
        utils.check_string_ids(controls)
        if not isinstance(perturbed, sc.AnnData):
            if isinstance(perturbed, list):
                if not all([isinstance(x, sc.AnnData) for x in perturbed]):
                    raise ValueError("Expected all perturbed datasets to be "
                                     "sc.AnnData objects.")
                if not all([utils.check_matching_genes(controls, x)\
                for x in perturbed]):
                    raise ValueError("Gene columns do not match between control"
                                     " and perturbed cells.")
                if not all([self.treatment_col in x.obs.columns\
                            for x in perturbed]):
                    raise ValueError("Expected {} column in perturbed data.".format(
                                     self.treatment_col))
                for each in perturbed:
                    utils.check_string_ids(each)
            else:
                raise ValueError("Unexpected input type for `perturbed`: "
                                 "{}. Expected list of sc.AnnData objects or "
                                 "a single sc.AnnData object".\
                                 format(type(perturbed)))
        else:
            if utils.check_matching_genes(controls, perturbed):
                raise ValueError("Gene columns do not match between control and"
                                    " perturbed cells.")
            if self.treatment_col not in perturbed.obs.columns:
                raise ValueError("Expected {} column in perturbed data.".format(
                                    self.treatment_col))
            utils.check_string_ids(perturbed)
        if self.reference == 'all':
            if verbose:
                print("Using all datasets as reference sets.")
            # logging.info('Using all datasets as references.')
            if isinstance(perturbed, list):
                reference = [controls.copy()] + [x.copy() for x in perturbed]
            else:
                reference = [controls.copy()] + [perturbed.copy()]
        else:
            if verbose:
                print("Using control samples as reference.")
            # logging.info("Using control samples as references.")
            reference = [controls.copy()]
        # if distance function is from ncfs.distances, expects feature weights 
        # check to see if they were provided, otherwise set weights to 1 
        if self.neighbor_kws['metric'].__module__ == 'ncfs.distances':
            try:
                self.neighbor_kws['metric_kwds']['w']
            except KeyError:
                # check for dict to update previously specified values
                if isinstance(self.neighbor_kws['metric_kwds'], dict):
                    self.neighbor_kws['metric_kwds'].update(
                                            {'w': np.ones(controls.X.shape[1])})
                else:
                    self.neighbor_kws['metric_kwds'] = {'w': np.ones(controls.X.shape[1])}
                    
        


        if self.cluster_col is not None:
            for adata in reference:
                if self.cluster_col not in adata.obs.columns:
                    raise ValueError(f"`cluster_col` - {self.cluster_col} not found"
                                    " in control data.")
                if np.any([pd.isnull(x) for x in adata.obs[self.cluster_col]]):
                    raise ValueError("Expected labelled cells by passing "\
                                     "`cluster_col`={}. ".format(self.cluster_col) +
                                     "Received atleast one unannotated cell.")

        # fit scale function to scale features between 0 and 1
        scaler = preprocessing.MinMaxScaler()
        if self.reference == 'all':
            scaler.fit(reference[0].concatenate(reference[1:], join='outer').X)
        else:
            scaler.fit(controls.X)
        if self.cluster_col is None:
            self.__cluster_references(reference, verbose)
        if verbose:
            print("-" * 20 + " Starting NCFS Fitting " + "-" * 20)
        # copy control clusters over to control dataset 
        controls.obs[self.cluster_col] = reference[0].obs[self.cluster_col]
        # logging.info("-" * 20 + " Starting NCFS Fitting " + "-" * 20)
        model, weights = self.__ncfs_fit(reference, scaler, verbose)
        if self.reference == 'all':
            treatments = [each.obs[self.treatment_col].values[0]\
                          for each in reference]
        del reference
        # logging.info('Removing reference data sets. Combining across treatments.')
        # utils.log_system_usage()
        # combine treatment datasets into single anndata object
        combined = controls.concatenate(perturbed, join='outer')
        # scale combined matrix between 0 and 1 and apply learned weights
        # across gene expression matrix
        combined.X = model.transform(scaler.transform(combined.X))
        # save genes weights
        combined.var['ncfs.weights'] = model.coef_
        if self.reference == 'all':
            for i, value in enumerate(treatments):
                combined.var["{}.weights".format(value)] = weights[i, :]
        n_clusters = len(controls.obs[self.cluster_col].unique())
        combined.var['informative'] = combined.var['ncfs.weights'].apply(lambda x: x > self.weight_threshold)
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
        sc.tl.umap(combined)
        # grab connectivities of cells
        if parse_version(sc.__version__) < parse_version("1.5.0"):
            A = combined.uns['neighbors']['connectivities']
        else:
            A = combined.obsp['connectivities']
        g = utils.igraph_from_adjacency(A)
        # instantiate semi-supervised Louvain model
        try:
            resolution = self.sslouvain_kws['resolution_parameter']
        except KeyError:
            resolution = 1.0
        try:
            vertex = self.sslouvain_kws['partition_type']
        except KeyError:
            vertex = sslouvain.RBConfigurationVertexPartition
        if not isinstance(vertex, utils.ig.VertexClustering):
            vertex = sslouvain.RBConfigurationVertexPartition

        y_, mutables = utils.format_labels(combined.obs[self.cluster_col])
        if verbose:
            print("Running semi-supervised louvain community detection")
        # logging.info("Runing sslouvain")
        part = sslouvain.find_partition(g,
                                        vertex,
                                        # sslouvain.CPMVertexPartition,
                                        # sslouvain.RBERVertexPartition,
                                        initial_membership=y_,
                                        mutable_nodes=mutables,
                                        resolution_parameter=resolution)
        # store new cluster labels in cell metadata
        combined.obs['sslouvain'] = part.membership
        combined.obs['sslouvain'] = combined.obs['sslouvain'].astype('category')
        if verbose:
            print("ICAT complete.")
        # utils.close_log()
        return combined

    def __cluster_references(self, reference, verbose):
        if verbose:
            print("Clustering reference datasets...")
            # logging.info('Clustering reference datasets.')
        # utils.log_system_usage()
        for i, adata in enumerate(reference):
            if verbose:
                print("Clustering cells in reference {}".format(i + 1))
                # logging.info("Clustering cells in reference {}".format(i + 1))
            utils.log_system_usage()
            sc.pp.pca(adata, **self.pca_kws)
            sc.pp.neighbors(adata, **self.neighbor_kws)
            sc.tl.umap(adata)
            if self.clustering == 'louvain':
                sc.tl.louvain(adata, **self.cluster_kws)
                self.cluster_col = 'louvain'
            elif self.clustering == 'leiden':
                sc.tl.leiden(adata, **self.cluster_kws)
                self.cluster_col = 'leiden'
        # logging.info('Cells clustered')
        # utils.log_system_usage()

    def __ncfs_fit(self, reference, scaler, verbose):
        # no previous clustering provided, cluster using louvain or leiden
        weights = np.zeros((len(reference), reference[0].shape[1]))
        for i, adata in enumerate(reference):
            # fit gene weights using control dataset
            if verbose:
                start = time.time()
                print(f"Starting NCFS gradient ascent for dataset {i + 1}")
            X_train = adata.X
            y_train = adata.obs[self.cluster_col].values
            if self.sub_sample:
                splits = train_test_split(adata.X,
                                          adata.obs[self.cluster_col],
                                          train_size=self.train_size,
                                          stratify=adata.obs[self.cluster_col])
                X_train = splits[0]
                y_train = splits[2].values
            # logging.info(f"Starting NCFS gradient ascent for dataset {i + 1}")
            # utils.log_system_usage()
            model = ncfs.NCFS(**self.ncfs_kws)
            model.fit(scaler.transform(X_train),
                      np.array(y_train),
                      sample_weights='balanced')
            weights[i, :] = model.coef_
            # logging.info(msg)
            if verbose:
                msg = "Dataset {} NCFS complete: {} to convergence".format(i + 1, 
                              utils.ftime(time.time() - start))
                print(msg)
        model.coef_ = np.max(weights, axis=0)
        return model, weights



