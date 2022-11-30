import itertools

import apricot
import ncfs
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from sklearn import metrics, preprocessing, model_selection, neighbors
from scipy import sparse
import sslouvain
from tqdm import tqdm

from . import utils


def optimize_louvain(
    adata,
    min_neighbors=3,
    max_neighbors=50,
    neighbor_step=2,
    min_res=0.6,
    max_res=1.2,
    res_step=0.02,
    rep="X_pca",
    return_best=True,
    verbose=False,
    num_cores=1,
):
    """
    Find optimal resolution and neighbor parameters for clustering control cells.

    Performs a grid search over resolution and neighbors values to optimize
    the Calinski-Harabasz score. The workflow assumes all desired pre-processing
    prior to clustering has already taken place.

    Parameters
    ----------
    adata : sc.AnnData
        Single cell dataset of control cells to cluster.
    min_neighbors : int, optional
        Minimum number of neigbhors to use when constructing the neighbor graph,
        by default 3.
    max_neighbors : int, optional
        Maximum number of neigbhors to use when constructing the neighbor graph,
        by default 3. by default 50
    neighbor_step : int, optional
        Step size for spanning range between `min_neighbors` and `max_neighbors`,
        by default 2.
    min_res : float, optional
        Minimum resolution parameter to use during Louvain clustering, by
        default 0.6.
    max_res : int, optional
        Maximum resolution parameter to use during Louvain clustering, by
        default 1.2.
    res_step : float, optional
        Step size for spanning range between `min_res` and `max_res`, by default
        0.02.
    rep : str, optional
        Data representation to use for neighbor graph construction. Default is
        'X_pca' and the principal components will be used. Values can either be
        'X' for processed gene expression matrix, or any data stored in
        `adata.obsm`.
    return_best : bool, optional
        Whether to only return the best n_neighbor and resolution pair for
        louvain clustering, by default True. Otherwise return all results as a
        pandas DataFrame
    verbose : bool, optional
        Whether to verbosely print progress, by default False.
    num_cores : int, optional
        Number of cores to use during grid search, by default 1.

    Returns
    -------
    pd.Dataframe | np.array
        Data frame of scores over tested parameter space, OR an
        (n_neighbor, resolution) array containing best performing parameter pair
        if `return_best == True`.
    """
    res = np.arange(min_res, max_res + res_step, res_step)
    n_neighbors = range(min_neighbors, max_neighbors + neighbor_step, neighbor_step)
    test = itertools.product(n_neighbors, res)

    if rep == "X_pca" and rep not in adata.obsm.keys():
        print(
            "WARNING: pca not run prior to louvain optimization. "
            "Running now with default parameters."
        )
        sc.pp.pca(adata)
    elif rep != "X" and rep not in adata.obsm.keys():
        raise ValueError("Data representation {rep} not found in adata.obsm.")

    def cluster_cells(p):
        n, r = p
        sc.pp.neighbors(adata, n_neighbors=n, use_rep=rep)
        sc.tl.louvain(
            adata,
            resolution=r,
            key_added="Cluster",
        )
        if adata.obs.Cluster.nunique() > 1:
            cal_hara = metrics.calinski_harabasz_score(
                adata.obsm[rep], adata.obs.Cluster
            )
        else:
            cal_hara = np.nan
        return n, r, cal_hara

    results = Parallel(n_jobs=num_cores, verbose=verbose)(
        delayed(cluster_cells)(p) for p in test
    )
    df = pd.DataFrame(
        results, columns=["n_neighbor", "resolution", "calinski_harabasz"], dtype=float
    )
    df.n_neighbor = df.n_neighbor.astype(int)

    if return_best:
        return df.loc[
            df.calinski_harabasz.idxmax(), ["n_neighbor", "resolution"]
        ].values
    return df


def optimize_ncfs(
    adata,
    y,
    n_neighbors=5,
    n_splits=3,
    sigma_vals=[0.5, 1, 1.5, 2, 2.5, 3],
    reg_vals=[0.25, 0.5, 1, 1.5, 2, 2.5, 3],
    max_cells=750,
    return_best=True,
    verbose=False,
):
    """Optimize NCFS regularization parameters.

    Optimize performance using k-fold cross validation for different values of
    NCFS `sigma` and `reg` parameters over a grid. Performance is evaluated
    using a KNN classifier with learned NCFS weights applied to highly variable
    genes in the provided dataset. Performance is measured using the
    Matthew's Correlation Co-efficient.

    Parameters
    ----------
    adata : sc.AnnData
        Dataset to assess NCFS parameters over.
    y : numpy.ndarray
        Class labels for cells in `adata`
    n_neighbors : int, optional
        Number of neighbors to use in knn classifier for evaluating performance.
    n_splits : int, optional
        Number of splits to use during k-fold cross validation. By default 3.
    sigma_vals : list[float], optional
        Kernel width hyperparameter values to test. By default, test 0.5, 1,
        1.5, 2, 2.5, 3.
    reg_vales : list[float], optional
        Regularization parameter values to test. By Default, 0.25, 0.5, 1, 1.5,
        2, 2.5, 3.
    max_cells : int, optional
        Maximum number of cells to train NCFS model on. Default is 750, and data
        sets with more cells will be selected down using submodular
        optimization.
    return_best : bool, optional
        Whether to only return the best sigma and reg pair for NCFS feature
        weighting, by default True. Otherwise return summarized results for each
        tested pair of values as a DataFrame.
    verbose : bool, optional
        Whether to verbosely print progress, by default False.

    Returns
    -------
    pd.Dataframe | np.array
        Data frame of scores over tested parameter space. Scores are summarized
        over each data split. Alternatively, a (sigma, reg) array containing
        best performing parameter pair if `return_best == True`.
    """
    adata = adata[:, adata.var.highly_variable].copy()
    if sparse.issparse(adata.X):
        print(
            "NCFS does not currently support sparse matrices, "
            "converting to dense."
        )
        adata.X = adata.X.todense()
    skf = model_selection.StratifiedKFold(n_splits=n_splits)
    scaler = preprocessing.MinMaxScaler().fit(adata.X)
    X = scaler.transform(adata.X)
    results = pd.DataFrame(
        itertools.product(sigma_vals, reg_vals, range(n_splits)),
        columns=["sigma", "reg", "split"],
    )
    results["mcc"] = np.nan
    # for sigma, reg in itertool.product(sigmavals, reg_vals)
    for i, (train_index, test_index) in enumerate(skf.split(X, y)):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # select cells to provide to NCFS using submodular optimization
        if train_index.size > max_cells:
            print(
                f"Train split > {max_cells} cells, performing submodular optimization "
                "to select representative cells"
            )
            X_train, y_train = utils.submodular_select(
                adata[train_index, :].copy(),
                y_train,
                by_cluster=False,
                stratified=True,
                train_size=max_cells / train_index.size,
                metric="euclidean",
                selector=apricot.FacilityLocationSelection,
                verbose=verbose,
            )
            X_train = scaler.transform(X_train)
        iterator = results.query(f"split=={i}").index
        if verbose:
            iterator = tqdm(iterator)
        for row in iterator:
            model = ncfs.NCFS(
                sigma=results.at[row, "sigma"], reg=results.at[row, "reg"]
            )
            model.fit(X_train, y_train)
            knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors).fit(
                model.transform(X_train), y_train
            )
            y_pred = knn.predict(model.transform(X_test))
            results.loc[row, "mcc"] = metrics.matthews_corrcoef(y_test, y_pred)
    summarized = (
        results.groupby(["sigma", "reg"]).median().drop(columns="split").reset_index()
    )
    if return_best:
        summarized.loc[summarized.mcc.idxmax(), ["sigma", "reg"]].values
    return summarized


def optimize_sslouvain(
    adata,
    control_clusters,
    reg,
    sigma,
    max_cells=750,
    min_neighbors=3,
    max_neighbors=50,
    neighbor_step=2,
    min_res=0.6,
    max_res=1.2,
    res_step=0.02,
    return_best=True,
    verbose=False,
    num_cores=1,
):
    """Optimize `n_neighbor` and `resolution` parameters for sslouvain.

    Performs a grid search over resolution and neighbor values to optimize
    the Calinski-Harabasz score. The workflow assumes all desired pre-processing
    prior to clustering has already taken place.


    Parameters
    ----------
    adata : sc.AnnData
        Single cell dataset of all cells to cluster. Should be N x M where N is
        the total number of cells in the dataset, and include control and
        treated cells.
    control_clusters : pd.Series
        Cluster labels for control cells  in `adata`. Should be N x 1 where N is
        the total number of cells in the dataset. Non-control cells should have
        null values (e.g. np.nan, None, etc.).
    reg : float
        Regulation parameter for NCFS. Can be found using `optimize_ncfs()`
    sigma : float
        Kernel width parameter for NCFS. Can be found using `optimize_ncfs()`
    min_neighbors : int, optional
        Minimum number of neigbhors to use when constructing the neighbor graph,
        by default 3.
    max_neighbors : int, optional
        Maximum number of neigbhors to use when constructing the neighbor graph,
        by default 3. by default 50
    neighbor_step : int, optional
        Step size for spanning range between `min_neighbors` and `max_neighbors`,
        by default 2.
    min_res : float, optional
        Minimum resolution parameter to use during Louvain clustering, by
        default 0.6.
    max_res : int, optional
        Maximum resolution parameter to use during Louvain clustering, by
        default 1.2.
    res_step : float, optional
        Step size for spanning range between `min_res` and `max_res`, by default
        0.02.
    rep : str, optional
        Data representation to use for neighbor graph construction. Default is
        'X_pca' and the principal components will be used. Values can either be
        'X' for processed gene expression matrix, or any data stored in
        `adata.obsm`.
    return_best : bool, optional
        Whether to only return the best n_neighbor and resolution pair for
        louvain clustering, by default True. Otherwise return all results as a
        pandas DataFrame
    verbose : bool, optional
        Whether to verbosely print progress, by default False.
    num_cores : int, optional
        Number of cores to use during grid search, by default 1.

    Returns
    -------
    _type_
        _description_
    """
    adata = adata[:, adata.var.highly_variable].copy()
    if sparse.issparse(adata.X):
        print(
            "NCFS does not currently support sparse matrices, ",
            "converting to dense."
        )
        adata.X = adata.X.todense()
    adata.X = preprocessing.MinMaxScaler().fit_transform(adata.X)
    control_index = ~control_clusters.isna()
    # subselect data
    if sum(control_index) > max_cells:
        X_train, y_train = utils.submodular_select(
            adata[control_index, :].copy(),
            control_clusters[control_index],
            by_cluster=False,
            stratified=True,
            train_size=max_cells / control_index.sum(),
            metric="euclidean",
            selector=apricot.FacilityLocationSelection,
            verbose=verbose,
        )
    else:
        X_train = adata[control_index, :].X
        y_train = control_clusters[control_index]
    if verbose:
        print("Fitting NCFS model...")
    model = ncfs.NCFS(
        reg=reg, sigma=sigma
    )
    model.fit(X_train, y_train)
    # multiply gene values by NCFS feature weights
    adata.X *= model.coef_
    # get set of values to test
    res = np.arange(min_res, max_res + res_step, res_step)
    n_neighbors = range(min_neighbors, max_neighbors + neighbor_step, neighbor_step)
    to_test = itertools.product(n_neighbors, res)

    def cluster_cells(p):
        n, r = p
        sc.pp.neighbors(adata, n_neighbors=n, use_rep='X')
        A = utils.get_neighbors(adata, "connectivities")
        g = utils.igraph_from_adjacency(A)
        vertex = sslouvain.RBConfigurationVertexPartition
        y_, mutables = utils.format_labels(control_clusters)
        part = sslouvain.find_partition(
            g,
            vertex,
            initial_membership=y_,
            mutable_nodes=mutables,
            resolution_parameter=r,
        )
        if len(set(part.membership)) > 1:
            cal_hara = metrics.calinski_harabasz_score(
                adata.X, part.membership
            )
        else:
            cal_hara = np.nan
        return n, r, cal_hara

    results = Parallel(n_jobs=num_cores, verbose=verbose)(
        delayed(cluster_cells)(p) for p in to_test
    )
    df = pd.DataFrame(
        results, columns=["n_neighbor", "resolution", "calinski_harabasz"], dtype=float
    )
    df.n_neighbor = df.n_neighbor.astype(int)

    if return_best:
        return df.loc[
            df.calinski_harabasz.idxmax(), ["n_neighbor", "resolution"]
        ].values
    return df