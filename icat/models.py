"""
Module to identify cell-types across treatments.

@author: Dakota Y. Hawkins
@contact: dyh0110@bu.edu

Example:

.. code-block::python
    from icat import simulate
    from icat import models
    import scanpy as sc
    import numpy as np
    data_model = simulate.SingleCellDataset(
        populations=2,
        genes=1000,
        dispersion=np.random.choice([1, 2, 3], 1000)
    )
    controls = data_model.simulate()
    controls.obs['treatment'] = 'control'
    perturbed = simulate.perturb(controls)
    perturbed.obs['treatment'] = 'perturbed'
    adata = controls.concatenate([perturbed])
    sc.pp.log1p(adata)
    model = models.icat(
        ctrl_value="control",
        ncfs_kws={'reg': 1, 'sigma': 3},
        neighbor_kws={'n_neighbors': 15}, 
        cluster_kws={'resolution': 0.75},
    )
    out = model.cluster(adata, adata.obs['treatment'])
    print(out.obs['sslouvain'].unique())
"""

import inspect
import warnings
import logging
import time

import numpy as np
import pandas as pd
from sklearn import preprocessing

from sklearn.neighbors import KNeighborsTransformer
import apricot
import scanpy as sc

import sslouvain
import ncfs

from . import utils


class icat:
    """
    Model to identify clusters across treatments.

    Parameters
    ----------
    ctrl_value : str
        Status label corresponding to control group in treatment observation
        data.
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
    subsample : bool, optional
        Whether to perform NCFS training on a sample of the data. Default is
        False, and training will be performed on the entire dataset.
    train_size : float, optional.
        Proportion of data to train on if `subsample==True`. Values should be
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
    neighbor_kws : dict, optional
        Keyword arguments for identifying neighbors. By default None, and
        default parameters are used. See `scanpy.pp.neighbors` for more
        information.
    sslouvain_kws : dict, optional
        Keyword arguments for `sslouvain`. By default None, and default
        parameters will be used. See `sslouvain.find_partition` for more
        information.
    pca_kws : dict, optional
        Keyword arguments for performing principle component analysis. By default
        None, and default parameters will be used. See `sslouvain.tl.pca` for
        more information.
    subsample_kws : dict, optional
        Keyword arguments for performing submodular optimization. By default
        None, and default parameters will be used. See `icat.select_cells` for
        more information.
    """

    def __init__(
        self,
        ctrl_value,
        clustering="louvain",
        reference="all",
        subsample=False,
        train_size=0.1,
        ncfs_kws=None,
        cluster_kws=None,
        cluster_col=None,
        neighbor_kws=None,
        sslouvain_kws=None,
        pca_kws=None,
        subsample_kws=None,
    ):
        self.ctrl_value = ctrl_value
        self.clustering = clustering
        self.reference = reference
        self.subsample = subsample
        self.train_size = train_size
        self.ncfs_kws = ncfs_kws
        self.cluster_kws = cluster_kws
        self.cluster_col = cluster_col
        self.neighbor_kws = neighbor_kws
        self.sslouvain_kws = sslouvain_kws
        self.pca_kws = pca_kws
        self.subsample_kws = subsample_kws
        # utils.set_log()

    @property
    def ctrl_value(self):
        return self._ctrl_value

    @ctrl_value.setter
    def ctrl_value(self, value):
        if not isinstance(value, str):
            raise ValueError("Expected string for `ctrl_value` parameter")
        self._ctrl_value = value

    @property
    def clustering(self):
        return self._clustering

    @clustering.setter
    def clustering(self, value):
        if not isinstance(value, str):
            raise ValueError(
                "Expected string for `clustering` parameter. "
                "Received {}".format(value)
            )
        value = value.lower()
        if value not in ["louvain", "leiden"]:
            raise ValueError("Unsupported method: {}".format(value))
        self._clustering = value

    @property
    def reference(self):
        """Reference data set. Acceptable values 'all' or 'controls'."""
        return self._reference

    @reference.setter
    def reference(self, value):
        if value not in ["all", "controls"]:
            raise ValueError("`reference` value should be either 'all' or 'controls'")
        self._reference = value

    @property
    def subsample(self):
        "Whether to perform NCFS training on a sample of the data."
        return self._subsample

    @subsample.setter
    def subsample(self, value):
        if not isinstance(value, bool):
            raise ValueError("Expected bool for `subsample`")
        self._subsample = value

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
        default_kws = utils.get_default_kwargs(ncfs.NCFS, ["self"])
        if value is not None:
            value = utils.check_kws(default_kws, value, "ncfs")
        else:
            value = default_kws
        self._ncfs_kws = value

    @property
    def neighbor_kws(self):
        """Keyword arguments for neighbor determination."""
        return self._neighbor_kws

    @neighbor_kws.setter
    def neighbor_kws(self, value):
        default_kws = utils.get_default_kwargs(sc.pp.neighbors, ["adata", "copy"])
        default_kws["metric"] = ncfs.distances.phi_s
        if value is not None:
            value = utils.check_kws(default_kws, value, "neighbor_kws")
        else:
            value = default_kws
        self._neighbor_kws = value

    @property
    def cluster_kws(self):
        """Keyword arguments to pass to clustering algorithm."""
        return self._cluster_kws

    @cluster_kws.setter
    def cluster_kws(self, value):
        default_kws = {
            "louvain": utils.get_default_kwargs(sc.tl.louvain, ["adata", "copy"])
        }
        if "leiden" in dir(sc.tl):
            default_kws["leiden"] = utils.get_default_kwargs(
                sc.tl.leiden, ["adata", "copy"]
            )
        else:
            if self.clustering == "leiden":
                raise ValueError(
                    "The version of Scanpy you are running does "
                    "not support Leiden clustering. If you wish "
                    "to use Leiden clustering, try upgrading your"
                    "version of Scanpy."
                )
        default_kws["resolution"] = 1.0
        if value is not None:
            value = utils.check_kws(
                default_kws[self.clustering], value, "cluster." + self.clustering
            )
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
            raise ValueError(
                "Expected integer index or string name for " + "`cluster_col`"
            )
        self._cluster_col = value

    @property
    def sslouvain_kws(self):
        """Keyword arguments for sslouvain."""
        return self._sslouvain_kws

    @sslouvain_kws.setter
    def sslouvain_kws(self, value):
        default_kws = utils.get_default_kwargs(sslouvain.find_partition, "")
        if value is not None:
            value = utils.check_kws(default_kws, value, "sslouvain_kws")
        else:
            value = default_kws
        self._sslouvain_kws = value

    @property
    def pca_kws(self):
        """Keyword arguments for PCA."""
        return self._pca_kws

    @pca_kws.setter
    def pca_kws(self, value):
        default_kws = utils.get_default_kwargs(sc.pp.pca, ["data", "copy"])
        if value is not None:
            value = utils.check_kws(default_kws, value, "pca_kws")
        else:
            value = default_kws
        self._pca_kws = value

    @property
    def subsample_kws(self):
        """Keyword arguments for `select_cells()`"""
        return self._subsample_kws

    @subsample_kws.setter
    def subsample_kws(self, value):
        default_kws = utils.get_default_kwargs(
            self.select_cells, ["adata", "label_col"]
        )
        if value is not None:
            if "selector" in value:
                if not any(
                    [
                        value["selector"] == x
                        for x in [
                            apricot.MaxCoverageSelection,
                            apricot.MixtureSelection,
                            apricot.SumRedundancySelection,
                            apricot.FacilityLocationSelection,
                            apricot.FeatureBasedSelection,
                            apricot.GraphCutSelection,
                            apricot.SaturatedCoverageSelection,
                        ]
                    ]
                ):
                    raise ValueError(
                        "Expected `apricot` submodular "
                        "optimization method for `selector`."
                    )
            value = utils.check_kws(default_kws, value, "subsample_kws")
        else:
            value = default_kws
        self._subsample_kws = value

    # @profile
    def cluster(self, adata, treatment, verbose=False):
        """
        Cluster cells in control and experimental conditions.

        Parameters
        ----------
        adata : sc.AnnData
            Annotated dataframe of cells.
        perturbed : pd.Series, numpy.ndarray, list-like
            Treatment labels for each cell in `adata`.
        verbose : boolean, optional
            Whether to print progress through clustering step.

        Returns
        -------
        sc.AnnData
            Annotated dataframe of out cells in NCFS space.
        """
        self.log_ = False
        if self.log_:
            utils.set_log()
            utils.log_system_usage("Instantiation")

        self.verbose_ = verbose
        if not isinstance(adata, sc.AnnData):
            raise ValueError("Expected AnnData object for `controls`.")
        if isinstance(treatment, pd.Series):
            treatment = treatment.values
        else:
            treatment = utils.check_np_castable(treatment, "treatment")
        if not self.ctrl_value in treatment:
            raise ValueError(
                f"Control value {self.ctrl_value} not found " "in treatment array."
            )
        # if distance function is from ncfs.distances, expects feature weights
        # check to see if they were provided, otherwise set weights to 1
        metric = self.neighbor_kws["metric"]
        if not isinstance(metric, str) and metric.__module__ == "ncfs.distances":
            try:
                self.neighbor_kws["metric_kwds"]["w"]
            except KeyError:
                # check for dict to update previously specified values
                if isinstance(self.neighbor_kws["metric_kwds"], dict):
                    self.neighbor_kws["metric_kwds"].update(
                        {"w": np.ones(adata.shape[1])}
                    )
                else:
                    self.neighbor_kws["metric_kwds"] = {"w": np.ones(adata.shape[1])}
        if self.log_:
            utils.log_system_usage("Starting NCFS transformation.")
        self.__learn_weights(adata, treatment)

        # n_clusters = len(controls.obs[self.cluster_col].unique())
        adata.var["informative"] = adata.var["ncfs_weights"].apply(lambda x: x > 1)
        informative = adata.var["informative"].sum()
        if self.verbose_:
            print(f"Found {informative} informative features.")

        # return adata
        # if informative < n_clusters:
        #     warnings.warn("Number of informative genes less "
        #                   "than the number of identified control clusters: "
        #                   f"informative genes: {informative}, "
        #                   f"number of clusters: {n_clusters}. Consider "
        #                   "increasing `sigma` or decreasing `reg` for better "
        #                   "performance.")
        # scikit-learn 0.22, umap==0.4.4
        self.neighbor_kws["use_rep"] = "X_icat"
        if self.log_:
            utils.log_system_usage("Before NCFS neighbors.")
        # A = KNeighborsTransformer(mode='connectivity',
        #                           n_neighbors=self.neighbor_kws['n_neighbors'])\
        #                          .fit_transform(adata.obsm['X_icat'])

        sc.pp.neighbors(adata, **self.neighbor_kws)
        # sc.pp.neighbors(adata, n_neighbors=self.neighbor_kws['n_neighbors'])
        # sc.pp.neighbors(adata)
        if self.log_:
            utils.log_system_usage("After NCFS neighbors.")
        # grab connectivities of cells
        A = utils.get_neighbors(adata, "connectivities")
        g = utils.igraph_from_adjacency(A)
        # instantiate semi-supervised Louvain model
        try:
            resolution = self.cluster_kws["resolution"]
            if resolution is None:
                resolution = 1
        except KeyError:
            resolution = 1.0
        try:
            vertex = self.sslouvain_kws["partition_type"]
        except KeyError:
            vertex = sslouvain.RBConfigurationVertexPartition
        if not isinstance(vertex, utils.ig.VertexClustering):
            vertex = sslouvain.RBConfigurationVertexPartition
        adata.obs.loc[
            [x != self.ctrl_value for x in treatment], self.cluster_col_
        ] = None
        y_, mutables = utils.format_labels(adata.obs[self.cluster_col_])
        if self.verbose_:
            print(f"{sum(mutables)} cells set as mutable out of {len(y_)}")
        # logging.info("Runing sslouvain")
        part = sslouvain.find_partition(
            g,
            vertex,
            initial_membership=y_,
            mutable_nodes=mutables,
            resolution_parameter=resolution,
        )
        # store new cluster labels in cell metadata
        adata.obs["sslouvain"] = part.membership
        adata.obs["sslouvain"] = adata.obs["sslouvain"].astype("category")
        if self.verbose_:
            print("ICAT complete.")
        # utils.close_log()
        return adata

    # @profile
    def __learn_weights(self, adata, treatment):
        reference = [utils.subset_cells(adata, treatment, self.ctrl_value)]
        scaler = preprocessing.MinMaxScaler()
        if self.reference == "all":
            if self.verbose_:
                print("Using all datasets as reference sets.")
            reference += [
                utils.subset_cells(adata, treatment, x, False)
                for x in np.unique(treatment)
                if x != self.ctrl_value
            ]
            scaler.fit(adata.X)
        else:
            if self.verbose_:
                print("Using control samples as reference.")
            scaler.fit(reference[0].X)
        reference = self.__cluster_references(reference)
        model, weights = self.__ncfs_fit(reference, scaler)
        if self.verbose_:
            print("-" * 20 + " Finished NCFS Fitting " + "-" * 20)

        # set default label to None, supply reference cluster labels
        adata.obs[self.cluster_col_] = None
        # match selected cells in split reference datasets to combined
        if self.subsample:
            adata.obs["selected"] = None
        for ref in reference:
            adata.obs.loc[ref.obs.index, self.cluster_col_] = ref.obs[self.cluster_col_]
            if self.subsample:
                adata.obs.loc[ref.obs.index, "selected"] = ref.obs["selected"]
        del reference
        adata.obsm["X_icat"] = model.transform(scaler.transform(adata.X))

        # copy control clusters over to control dataset

        # save genes weights
        adata.var["ncfs_weights"] = model.coef_
        if self.reference == "all":
            for i, value in enumerate(np.unique(treatment)):
                adata.var["{}_weights".format(value)] = weights[i, :]
        # free up memory by deleting unecessary objects
        del model
        try:
            del adata.uns["neighbors"]
        except (AttributeError, KeyError):
            pass
        try:
            del adata.obsp["neighbors"]
        except (AttributeError, KeyError):
            pass
        if self.log_:
            utils.log_system_usage("NCFS transform complete.")

    # @profile
    def __cluster_references(self, reference):
        # logging.info('Clustering reference datasets.')
        # utils.log_system_usage()
        for i, ref in enumerate(reference):
            if self.verbose_:
                print("Clustering cells in reference {}.".format(i + 1))
                # logging.info("Clustering cells in reference {}".format(i + 1))
            # utils.log_system_usage()
            sc.pp.pca(ref, **self.pca_kws)
            sc.pp.neighbors(ref, **self.neighbor_kws)
            if self.cluster_col is not None:
                print(f"Using provided cell labels in {self.cluster_col}.")
                if self.cluster_col not in ref.obs.columns:
                    raise ValueError(
                        f"Provided cluster column {self.cluster_col} "
                        "not found in observation data."
                    )
                ref.obs[self.cluster_col + "-copy"] = ref.obs[self.cluster_col].copy()
                self.cluster_col_ = self.cluster_col + "-copy"
            else:
                if self.clustering == "louvain":
                    sc.tl.louvain(ref, **self.cluster_kws)
                    self.cluster_col_ = "louvain"
                elif self.clustering == "leiden":
                    sc.tl.leiden(ref, **self.cluster_kws)
                    self.cluster_col_ = "leiden"
        if self.log_:
            utils.log_system_usage("Clustered reference datasets.")
        return reference
        # logging.info('Cells clustered')
        # utils.log_system_usage()

    # @profile
    def select_cells(
        self,
        adata,
        label_col,
        method="submodular",
        selector=apricot.FacilityLocationSelection,
        use_rep="X",
        by_cluster=False,
        stratified=True,
    ):
        """
        Select cells to train NCFS weights on.

        Parameters
        ----------
        adata : sc.AnnData
            Annotated dataframe containing cell labels and expression values.
        label_col : str
            Column in `adata.obs` containing labels to use during NCFS.
        method : str, optional
            Method to select training cells. Possible values are 'submodular'
            and 'random'. By default 'submodular', and sobmodular optimization
            will be used to find the "best" cells to train on.
        selector : apricot.Function, optional
            Sub-optimal optimizer function from `apricot`. Default is
            `apricot.MaxCoverageSelection`.
        by_cluster : bool, optional
            Whether to perform submodular optimization on a per cluster basis.
            Otherwise perform on the whole dataset with labels provided. By
            default true.
        stratified : bool, optional
            Whether to select cells in a stratified manner. By default False,
            and for `by_cluster` submodular optimization, an equal number of
            cells from each cluster will be chosen. Otherwise, a number
            proportional to the label prevalence will be chosen.

        Returns
        -------
        (X_train, y_train) : (np.ndarray, np.ndarray)
            A tuple of the data matrix X and observation labels y to be used
            during NCFS training.
        """

        y = adata.obs[label_col].values
        if method == "submodular":
            train_X, train_y = utils.submodular_select(
                adata,
                y,
                by_cluster,
                stratified,
                self.train_size,
                self.neighbor_kws["metric"],
                selector,
                self.neighbor_kws["metric_kwds"],
                use_rep=use_rep,
                verbose=self.verbose_,
            )
        elif method == "random":
            train_X, train_y = utils.random_select(
                adata,
                y,
                self.train_size,
                stratified,
                use_rep=use_rep,
                verbose=self.verbose_,
            )
        elif method == "centroid":
            train_X, train_y = utils.centroid_collapse(
                adata, y, use_rep=use_rep, verbose=self.verbose_
            )

        elif method == "centroid_knn":
            train_X, train_y = utils.centroid_neighbors(
                adata,
                y,
                self.train_size,
                self.neighbor_kws["metric"],
                stratified,
                self.neighbor_kws["metric_kwds"],
                use_rep=use_rep,
                verbose=self.verbose_,
            )
        else:
            raise NotImplementedError(
                f"No support for method {method}. "
                "Supported methods are 'submodular' and "
                "random'."
            )
        if self.log_:
            utils.log_system_usage("After cell selection.")
        return train_X, train_y

    # @profile
    def __ncfs_fit(self, reference, scaler):
        # no previous clustering provided, cluster using louvain or leiden
        weights = np.zeros((len(reference), reference[0].shape[1]))
        for i, ref in enumerate(reference):
            X_train = ref.X
            y_train = ref.obs[self.cluster_col_].values
            if self.subsample:
                X_train, y_train = self.select_cells(
                    ref, self.cluster_col_, **self.subsample_kws
                )
            if self.verbose_:
                print(
                    f"Training NCFS on {X_train.shape[0]} samples out of "
                    + f"{ref.shape[0]}"
                )
            X_train = scaler.transform(X_train)
            # fit gene weights using control dataset
            if self.verbose_:
                start = time.time()
                print(
                    f"Starting NCFS gradient ascent for dataset {i + 1}/{len(reference)}"
                )

            # logging.info(f"Starting NCFS gradient ascent for dataset {i + 1}")
            # utils.log_system_usage()
            y_train, __ = utils.format_labels(y_train)
            model = ncfs.NCFS(**self.ncfs_kws)
            model.fit(X_train, y_train, sample_weights="balanced")
            weights[i, :] = model.coef_
            # logging.info(msg)
            if self.verbose_:
                msg = "Dataset {} NCFS complete: {} to convergence".format(
                    i + 1, utils.ftime(time.time() - start)
                )
                print(msg)
        model.coef_ = np.max(weights, axis=0)
        if self.log_:
            utils.log_system_usage("After NCFS fitting.")
        return model, weights
