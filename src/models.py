
from scanpy import api as sc
import numpy as np
from sklearn import preprocessing, discriminant_analysis, neighbors
from ncfs_expanded import NCFS
import pandas as pd

from ssLouvain import ssLouvain

class icat():
    """
    Model to Identify Clusters Across Treatments. 
    """

    def __init__(self, method='ncfs', treatment_col='treatment',
                 method_kws=None, weight_threshold=None, neighbor_kws=None,
                 leiden_kws=None, sslouvain_kws=None):
        self.method = method
        self.treatment_col = treatment_col
        self.method_kws = method_kws
        self.weight_threshold = weight_threshold
        self.neighbor_kws = neighbor_kws
        self.leiden_kws = leiden_kws
        self.sslouvain_kws = sslouvain_kws
    
    @property
    def method(self):
        return self._method
    
    @method.setter
    def method(self, value):
        if not isinstance(value, str):
            raise ValueError("Expected string for `method` parameter." 
                             "Received {}".format(value))
        value = value.lower()
        if value not in ['ncfs', 'lda', 'qda']:
            raise ValueError("Unsupported method: {}".format(value))
        self._method = value

    @property
    def treatment_col(self):
        return self._treatment_col

    @treatment_col.setter
    def treatment_col(self, value):
        self._treatment_col = value

    @property
    def method_kws(self):
        return self._method_kws
    
    @method_kws.setter
    def method_kws(self, value):
        default_kws = {'ncfs': {'alpha': 0.01, 'metric': 'euclidean', 'reg': 3,
                                'sigma': 2, 'eta': 10e-6},
                       'lda': {'solver': 'eigen', 'shrinkage': 'auto',
                               'priors': None, 'n_components': None,
                               'store_covariance': False, 'tol': 0.0001},
                       'qda': {'solver': 'eigen', 'shrinkage': 'auto',
                               'priors': None, 'n_components': None,
                               'store_covariance': False, 'tol': 0.0001}}
        if value is not None:
            value = __check_kws(default_kws[self.method], value,
                                'method.' + self.method)
        else:
            value = default_kws[self.method]
        self._method_kws = value
    
    @property
    def weight_threshold(self):
        return self._weight_threshold
    
    @weight_threshold.setter
    def weight_threshold(self, value):
        if self.method != 'ncfs':
            print("Warning: `weight_threshold` only pertains to `ncfs` method. "
                  "Changing the value will have no affect on outcome.")
        if value is None:
            value = 1
        self._weight_threshold = value

    @property
    def neighbor_kws(self):
        return self._neighbor_kws

    @neighbor_kws.setter
    def neighbor_kws(self, value):
        default_kws = {'n_neighbors': 15, 'n_pcs': None, 'use_rep': None,
                       'knn': True, 'random_state': 0, 'method': 'umap',
                       'metric': 'euclidean', 'metric_kwds':{}}
        if value is not None:
            value = __check_kws(default_kws, value, 'neighbor_kws')
        else:
            value = default_kws
        self._neighbor_kws = value
        
    @property
    def leiden_kws(self):
        return self._leiden_kws
    
    @leiden_kws.setter
    def leiden_kws(self, value):
        default_kws = {'resolution': 1, 'random_state': 0,
                       'adjacency': None, 'directed': True, 'use_weights': True,
                       'n_iterations': -1, 'partition_type': None}
        if value is not None:
            value = __check_kws(default_kws, value, 'leiden_kws')
        else:
            value = default_kws
        self._leiden_kws = value

    @property
    def sslouvain_kws(self):
        return self._sslouvain_kws
    
    @sslouvain_kws.setter
    def sslouvain_kws(self, value):
        default_kws = {'immutable': True, 'precluster': False, 'seed': None}
        if value is not None:
            value = __check_kws(default_kws, value, 'sslouvain_kws')
        else:
            value = default_kws
        self._sslouvain_kws = value

    def cluster(self, controls, perturbed):
        if not isinstance(controls, sc.AnnData):
            raise ValueError("Expected AnnData objecto for `controls`.")
        if not isinstance(perturbed, sc.AnnData):
            if isinstance(perturbed, list):
                if not all([isinstance(x, sc.AnnData) for x in perturbed]):
                    raise ValueError("Expected all perturbed datasets to be "
                                    "sc.AnnData objects.")
                if not all([check_matching_genes(controls, x)\
                for x in perturbed]):
                    raise ValueError("Gene columns do not match between control"
                                     " and perturbed cells.")
                perturbed = sc.AnnData(
                                X=np.vstack([each.X for each in perturbed]),
                                obs=pd.concat([each.obs for each in perturbed]),
                                var=controls.var)
            else:
                raise ValueError("Unexpected input type for `perturbed`: "
                                "{}. Expected list of sc.AnnData objects or "
                                "a single sc.AnnData object".\
                                format(type(perturbed)))
        if check_matching_genes(controls, perturbed):
            raise ValueError("Gene columns do not match between control and"
                                " perturbed cells.")
        if self.treatment_col not in perturbed.obs.columns:
            raise ValueError("Expected {} column in perturbed data.".format(
                                self.treatment_col))
        # scale perturbed data using control data
        scaler = preprocessing.StandardScaler()
        scaler.fit(controls.X)
        perturbed = sc.AnnData(X=scaler.transform(perturbed.X),
                               obs=perturbed.obs, var=perturbed.var)
            
        # scale cells to 0 centered with unit variance
        sc.pp.scale(controls)
        sc.pp.neighbors(controls, **self._neighbor_kws)
        sc.tl.leiden(controls, **self._leiden_kws)

        if self.method == 'ncfs':
            model = NCFS.NCFS(**self.method_kws)
        elif self.method == 'lda':
            model = discriminant_analysis.LinearDiscriminantAnalysis(
                                          **self.method_kws)
        else:
            model = discriminant_analysis.QuadraticDiscriminantAnalysis(
                                          **self.method_kws)
        model.fit(controls.X, np.array(controls.obs['leiden'].values))
        X_ = model.transform(np.vstack((controls.X, perturbed.X)))
        if self.method == 'ncfs':
            selected = np.where(model.coef_ > self.weight_threshold)[0]
            X_ = X_[:, selected]
        A_ = neighbors.kneighbors_graph(X_,
                                  n_neighbors=self._neighbor_kws['n_neighbors'],
                                  mode='distance')
        ss_model = ssLouvain.ssLouvain(**self.sslouvain_kws)
        y_ = np.hstack([controls.obs['leiden'].values,
                        np.array([np.nan]*perturbed.shape[0])])
        ss_model.fit(A_, y_)
        out = sc.AnnData(X=np.vstack((controls.X, perturbed.X)),
                         obs=pd.concat([controls.obs, perturbed.obs]),
                         var=controls.var)
        out.obs['sslouvain'] = ss_model.labels_
        out.uns['cluster_dims'] = X_
        return out


def __check_kws(reference_dict, new_dict, name):
    if not isinstance(new_dict, dict):
        raise ValueError("Expected dictionary of keyword arguments for "
                         "`{}`. Received {}.".format(name, type(new_dict)))
    for key, item in new_dict.items():
        if key not in reference_dict.keys():
            raise ValueError("Unexpected keyword argument `{}` for "
                             "{} keywords.".format(key, name))
        new_dict[key] = item
    return new_dict


def check_matching_genes(ref, new):
    return set(ref.var.index.values).difference(new.var.index.values) == 0