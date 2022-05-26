# icat
Identifying Cell-states Across Treatments

ICAT is a tool developed to better identify cell states in scRNAseq experiments where perturbations or some other biologic heterogeneity is present, such as gene knock-outs.

The method works by first identifying a set of conrol-defined cell states by performing unsupervised clustering. These identified cell states are then fed into a sparse gene weighting algorithm, Neighborhood Component Feature Selection (NCFS), to highly weight the most predictive genes, while also removing variance from non-explanatory genes. We then transform the data matrix using this weight vector, and perform semi-supervised clustering such that the originally identified control labels remain constant, but cells from experimental conditions

## Instalation

ICAT is available on `pip` and can be installed using the following command:

`pip install icat`

## How to use

ICAT makes heavy use of the excellent `scanpy` library along with the associated `anndata` structure.

Assuming an `sc.AnnData` object, `adata`, ICAT can be run with following code:


```python
from icat import models
hvgs = adata[:, adata.var.highly_variable].copy()
# create `icat` object, set "control" as label defining control cells
icat_model = models.icat(ctrl_value="control", reference="control")
# pass expression data, as well as a vector containing treatment status for each cell in your dataset
out = icat_model.cluster(hvgs, adata.obs["treatment"], verbose=True)
# final cluster labels in newly added "sslouvain" column
print(out.obs.sslouvain)
# highly weighted genes used to cluster cells (NCFS weight > 1)
out.obs.query("informative')

```


