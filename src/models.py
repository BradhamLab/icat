import inspect

import numpy as np
import pandas as pd
from sklearn import preprocessing, neighbors
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA

from ncfs_expanded import NCFS
from scanpy import api as sc
from ssLouvain import ssLouvain

from icat.src import utils


class icat():
    """
    Model to Identify Clusters Across Treatments. 
    """

    def __init__(self, method='ncfs', clustering='louvain',
                 treatment_col='treatment', method_kws=None, cluster_kws=None,
                 cluster_col=None, weight_threshold=None, neighbor_kws=None,
                 sslouvain_kws=None, pca_kws=None, use_X='X'):
        self.method = method
        self.clustering = clustering
        self.treatment_col = treatment_col
        self.method_kws = method_kws
        self.cluster_kws = cluster_kws
        self.cluster_col = cluster_col
        self.weight_threshold = weight_threshold
        self.neighbor_kws = neighbor_kws
        self.sslouvain_kws = sslouvain_kws
        self.pca_kws = pca_kws
        self.use_X = use_X
    
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
        default_kws = {'ncfs': utils.get_default_kwargs(NCFS.NCFS, ['self']),
                       'lda': utils.get_default_kwargs(LDA, ['self']),
                       'qda': utils.get_default_kwargs(QDA, ['self'])}
        if value is not None:
            value = utils.check_kws(default_kws[self.method], value,
                                'method.' + self.method)
        else:
            value = default_kws[self.method]
        self._method_kws = value
    
    @property
    def weight_threshold(self):
        return self._weight_threshold
    
    @weight_threshold.setter
    def weight_threshold(self, value):
        if self.method != 'ncfs' and value is not None:
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
        default_kws = utils.get_default_kwargs(sc.pp.neighbors, ['adata', 'copy'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'neighbor_kws')
        else:
            value = default_kws
        self._neighbor_kws = value
        
    @property
    def cluster_kws(self):
        return self._cluster_kws
    
    @cluster_kws.setter
    def cluster_kws(self, value):
        default_kws = {'louvain': utils.get_default_kwargs(sc.tl.louvain,
                                                     ['adata', 'copy']),
                       'leiden': utils.get_default_kwargs(sc.tl.leiden,
                                                    ['adata', 'copy'])}
        if value is not None:
            value = utils.check_kws(default_kws[self.clustering],
                                    value, 'cluster.' + self.clustering)
        else:
            value = default_kws[self.clustering]
        self._cluster_kws = value

    @property
    def sslouvain_kws(self):
        return self._sslouvain_kws
    
    @sslouvain_kws.setter
    def sslouvain_kws(self, value):
        default_kws = utils.get_default_kwargs(ssLouvain.ssLouvain, 'self')
        if value is not None:
            value = utils.check_kws(default_kws, value, 'sslouvain_kws')
        else:
            value = default_kws
        self._sslouvain_kws = value

    @property
    def pca_kws(self):
        return self._pca_kws

    @pca_kws.setter
    def pca_kws(self, value):
        default_kws = utils.get_default_kwargs(sc.pp.pca, ['data', 'copy'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'pca_kws')
        else:
            value = default_kws
        self._pca_kws = value

    @property
    def use_X(self):
        return self._use_X

    @use_X.setter
    def use_X(self, value):
        if not isinstance(value, str):
            raise ValueError("Expected string for `use_X`. Got {}".\
                             format(type(value)))
        if value not in ['X', 'pca']:
            raise ValueError("Expected `X` or `pca`. Got {}".format(value))
        self._use_X = value

    def cluster(self, controls, perturbed):
        if not isinstance(controls, sc.AnnData):
            raise ValueError("Expected AnnData objecto for `controls`.")
        if not isinstance(perturbed, sc.AnnData):
            if isinstance(perturbed, list):
                if not all([isinstance(x, sc.AnnData) for x in perturbed]):
                    raise ValueError("Expected all perturbed datasets to be "
                                     "sc.AnnData objects.")
                if not all([utils.check_matching_genes(controls, x)\
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
        if utils.check_matching_genes(controls, perturbed):
            raise ValueError("Gene columns do not match between control and"
                                " perturbed cells.")
        if self.treatment_col not in perturbed.obs.columns:
            raise ValueError("Expected {} column in perturbed data.".format(
                                self.treatment_col))
        # scale perturbed data using control data
        scaler = preprocessing.StandardScaler()
        scaler.fit(controls.X)
        sc.pp.pca(controls, **self.pca_kws)
        sc.pp.neighbors(controls, **self.neighbor_kws)
        sc.tl.umap(controls, min_dist=0.0)
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
        # scale cells to 0 centered with unit variance
        fit_X = scaler.transform(controls.X).astype(np.float64)
        perturb_X = scaler.transform(perturbed.X).astype(np.float64)
        if self.use_X == 'pca':
            pca_model = PCA(n_components=self.pca_kws['n_comps'],
                            svd_solver=self.pca_kws['svd_solver'],
                            random_state=self.pca_kws['random_state'])
            fit_X = pca_model.fit_transform(controls.X)
            perturb_X = pca_model.transform(perturb_X)
        
        model.fit(fit_X, np.array(controls.obs[cluster_col].values))
        X_ = model.transform(np.vstack((fit_X, perturb_X)))
        if self.method == 'ncfs':
            selected = np.where(model.coef_ > self.weight_threshold)[0]
            if len(selected) == 0:
                print('WARNING: No feature weights met threshold criteria. '
                      'All genes will be used. Try lowering threshold value for'
                      ' future runs.')
                selected = np.arange(len(model.coef_))
            X_ = X_[:, selected]
            var_ = controls.var.iloc[selected, :]
        elif self.method == 'lda':
            var_ = pd.DataFrame(['LDA.{}'.format(i + 1)\
                               for i in range(X_.shape[1])],
                               columns=['Dimension'])
        else:
            var_ = pd.DataFrame(['QDA.{}'.format(i + 1)\
                    for i in range(X_.shape[1])],
                    columns=['Dimension'])
        combined = sc.AnnData(X=X_,
                              obs=pd.concat([controls.obs, perturbed.obs],
                                            axis=0,
                                            sort=False).reset_index(drop=True),
                              var=var_)
        sc.pp.neighbors(combined, **self.neighbor_kws)
        sc.tl.umap(combined, min_dist=0.0)
        # if combined.shape[1] > 10:
        #     A_ = combined.uns['neighbors']['connectivities']
        # else:
        #     A_ = neighbors.kneighbors_graph(X_,
        #                                     self.neighbor_kws['n_neighbors'],
        #                                     mode='distance')
        A_ = combined.uns['neighbors']['connectivities']
        ss_model = ssLouvain.ssLouvain(**self.sslouvain_kws)
        y_ = np.hstack([controls.obs[cluster_col].values,
                        np.array([np.nan]*perturbed.shape[0])])
        ss_model.fit(A_, y_)

        combined.obs['sslouvain'] = ss_model.labels_
        return combined


if __name__ == '__main__':
    import sys
    sys.path.append('src/')
    import simulate
    data_model = simulate.SingleCellDataset()
    controls = data_model.simulate()
    perturbed = simulate.perturb(controls)
    perturbed.obs['treatment'] = 'ayo'
    model = icat(method='lda', method_kws={'shrinkage': 0.75, 'solver':'eigen'},
                 neighbor_kws={'n_neighbors': 100},
                 sslouvain_kws={'immutable': True, 'precluster': True},
                 cluster_kws={'resolution': 1.25})
    out = model.cluster(controls, perturbed)
