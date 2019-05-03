import inspect

import numpy as np
import pandas as pd
from sklearn import neighbors, preprocessing
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA

from ncfs_expanded import NCFS
from scanpy import api as sc
from ssLouvain import ssLouvain


class icat():
    """
    Model to Identify Clusters Across Treatments. 
    """

    def __init__(self, method='ncfs', clustering='louvain',
                 treatment_col='treatment', method_kws=None, cluster_kws=None,
                 cluster_col=None, weight_threshold=None,
                 neighbor_kws=None, sslouvain_kws=None):
        self.method = method
        self.clustering = clustering
        self.treatment_col = treatment_col
        self.method_kws = method_kws
        self.cluster_kws = cluster_kws
        self.cluster_col = cluster_col
        self.weight_threshold = weight_threshold
        self.neighbor_kws = neighbor_kws
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
        return self._treatment_col

    @treatment_col.setter
    def treatment_col(self, value):
        self._treatment_col = value

    @property
    def method_kws(self):
        return self._method_kws
    
    @method_kws.setter
    def method_kws(self, value):
        default_kws = {'ncfs': get_default_kwargs(NCFS.NCFS, ['self']),
                       'lda': get_default_kwargs(LDA, ['self']),
                       'qda': get_default_kwargs(QDA, ['self'])}
        if value is not None:
            value = _check_kws(default_kws[self.method], value,
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
        default_kws = get_default_kwargs(sc.pp.neighbors, ['adata', 'copy'])
        if value is not None:
            value = _check_kws(default_kws, value, 'neighbor_kws')
        else:
            value = default_kws
        self._neighbor_kws = value
        
    @property
    def cluster_kws(self):
        return self._cluster_kws
    
    @cluster_kws.setter
    def cluster_kws(self, value):
        default_kws = {'louvain': get_default_kwargs(sc.tl.louvain,
                                                     ['adata', 'copy']),
                       'leiden': get_default_kwargs(sc.tl.leiden,
                                                    ['adata', 'copy'])}
        if value is not None:
            value = _check_kws(default_kws, value, 'cluster.' + self.cluster)
        else:
            value = default_kws[self.clustering]
        self._cluster_kws = value

    @property
    def sslouvain_kws(self):
        return self._sslouvain_kws
    
    @sslouvain_kws.setter
    def sslouvain_kws(self, value):
        default_kws = get_default_kwargs(ssLouvain.ssLouvain, 'self')
        if value is not None:
            value = _check_kws(default_kws, value, 'sslouvain_kws')
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
        sc.pp.neighbors(controls, **self.neighbor_kws)
        if self.clustering == 'louvain':
            sc.tl.louvain(controls, **self.cluster_kws)
            cluster_col = 'louvain'
        elif self.clustering == 'leiden':
            sc.tl.leiden(controls, **self.cluster_kws)
            cluster_col = 'leiden'
        elif self.cluster_col is not None:
            cluster_col = self.cluster_col

        if self.method == 'ncfs':
            model = NCFS.NCFS(**self.method_kws)
        elif self.method == 'lda':
            model = LDA(**self.method_kws)
        else:
            model = QDA(**self.method_kws)
        model.fit(np.array(controls.X, dtype=np.float64),
                  np.array(controls.obs[cluster_col].values))
        X_ = model.transform(np.vstack((controls.X, perturbed.X)))
        if self.method == 'ncfs':
            selected = np.where(model.coef_ > self.weight_threshold)[0]
            if len(selected) == 0:
                print('WARNING: No feature weights met threshold criteria. '
                      'All genes will be used. Try lowering threshold value for'
                      ' future runs.')
                selected = np.arange(len(model.coef_))
            X_ = X_[:, selected]
        A_ = neighbors.kneighbors_graph(X_,
                                  n_neighbors=self._neighbor_kws['n_neighbors'],
                                  mode='distance').toarray()
        ss_model = ssLouvain.ssLouvain(**self.sslouvain_kws)
        y_ = np.hstack([controls.obs[cluster_col].values,
                        np.array([np.nan]*perturbed.shape[0])])
        ss_model.fit(A_, y_)
        out = sc.AnnData(X=np.vstack((controls.X, perturbed.X)),
                         obs=pd.concat([controls.obs, perturbed.obs]),
                         var=controls.var)
        out.obs['sslouvain'] = ss_model.labels_
        out.uns['cluster_dims'] = X_
        return out


def _check_kws(reference_dict, new_dict, name):
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


def get_default_kwargs(func, ignore_params):
    params = inspect.signature(func).parameters
    kwargs = {x:params[x].default for x in params if x not in ignore_params}
    return kwargs


if __name__ == '__main__':
    import sys
    sys.path.append('src/')
    import simulate
    data_model = simulate.SingleCellDataset()
    controls = data_model.simulate()
    perturbed = simulate.perturb(controls)
    perturbed.obs['treatment'] = 'ayo'
    model = icat(method_kws={'kernel': 'gaussian', 'metric': 'euclidean'},
                 neighbor_kws={'n_neighbors': 100})
    out = model.cluster(controls, perturbed)
