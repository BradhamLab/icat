# icat
Identifying Cell-states Across Treatments

ICAT is a tool developed to better identify cell states in scRNAseq experiments where perturbations or some other biologic heterogeneity is present, such as gene knock-outs.

The method works by first identifying a set of conrol-defined cell states by performing unsupervised clustering. These identified cell states are then fed into a sparse gene weighting algorithm, Neighborhood Component Feature Selection (NCFS), to highly weight the most predictive genes, while also removing variance from non-explanatory genes. We then transform the data matrix using this weight vector, and perform semi-supervised clustering such that the originally identified control labels remain constant, but cells from experimental conditions are free to cluster with any other cells regardless of treatment status.

## Installation

ICAT can be installed on linux machines using `pip` with the following command:

`pip install icat-sc`

## Pre-print
To learn more about the algorithm, and how it compares to other methods, see our pre-print on [BioArxiv](https://www.biorxiv.org/content/10.1101/2022.05.26.493603v2)

## How to use

ICAT makes heavy use of the excellent `scanpy` library along with the associated `AnnData` data structure.

An example code block walks through running `icat` on a simulated dataset. The 
final clustering is stored in the `sslouvain` column of the returned `AnnData`
object.

```python
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
```
**visualizing dataset**
![](docs/images/raw_input.png)
```python
    # specify model parameters -- see documentation for more information
    model = models.icat(
        ctrl_value="control",
        ncfs_kws={'reg': 1, 'sigma': 3},
        neighbor_kws={'n_neighbors': 15}, 
        cluster_kws={'resolution': 0.75},
    )
    # cluster cells by providing treatment information
    out = model.cluster(adata, adata.obs['treatment'])
    print(out.obs['sslouvain'].unique())
```
**visualizing results**

While ICAT does not automatically compute UMAP, tsne, or other reduced dimension
visualizations during clustering, it is possible to pass the upweighted count matrix
(found in `adata.obsm["X_icat"]`) to these algorithms. In the case of UMAP, the
returned `adata` object already has neighbors defined in this upweighted space, so 
calculating a new UMAP is simple:

```
sc.tl.umap(out)
sc.pl.umap(out, colors=['sslouvain', 'Population'])
```

![](docs/images/icat_output.png)

## Hyper Parameter Optimization
For working with your own data, we recommend finding appropriate Louvain and NCFS hyper
parameters **prior** to running the complete ICAT workflow. All hyper parameters used in the
original pre-print can be found as supplemental tables.

We have also provided grid search functions to find the "best" `n_neighbor` and `resolution`
parameters for Louvain and Semi-supervised Louvain clustering steps, as well as a function
to find the "best" kernel width (`sigma`) and regularization parameters (`reg)`. 

```python
from icat import optimize
sc.pp.pca(controls)

# Find the "best" `n_neighbor` and `resolution` parameter for clustering control cells by
# optimizing the Calinski-Harabasz Index over a grid of `n` and `r` values
louvain_n, louvain_r = optimize.optimize_louvain(
    controls,
    min_neighbors=3,
    max_neighbors=50,
    neighbor_step=2,
    min_res=0.3,
    max_res=1.2,
    res_step=0.02,
)

# cluster control cells with "best" values
sc.pp.neighbors(controls, n_neighbors=louvain_n)
sc.tl.louvain(controls, resolution=louvain_r)

# find "best" `sigma` and `reg` NCFS values by measuring the MCC in a
# weighted KNN over k-fold cross validation
sigma, reg =  optimize.optimize_ncfs(
    controls,
    controls.obs.louvain,
    n_neighbors=5,
    n_splits=3,
    sigma_vals=[0.5, 1, 1.5, 2, 2.5, 3],
    reg_vals=[0.25, 0.5, 1, 1.5, 2, 2.5, 3],
)
```

By default, ICAT uses the same `n_neighbors` and `resolution` parameters during
semi-supervised clustering as it does during control clustering. In practice,
this leads to good results (see paper). However, if users _would_ like to
optimize these parameters separately, we've include the below function:

```python
# optimize louvain parameters for semi-supervised clustering in NCFS
# space -- include complete dataset now 
adata.obs.loc[controls.obs.index, "control_clusters"] = controls.obs.louvain
sslouvain_n, sslouvain_r = optimize_sslouvain(
    adata,
    adata.obs.control_clusters,
    reg,
    sigma,
    max_cells=750,
    min_neighbors=3,
    max_neighbors=50,
    neighbor_step=2,
    min_res=0.3,
    max_res=1.2,
    res_step=0.02,
)
```
