import inspect

import numpy as np
import pandas as pd
from sklearn import preprocessing, neighbors
from sklearn.exceptions import NotFittedError
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
from scipy import stats

from ncfs_expanded import NCFS
from scanpy import api as sc
from sslouvain import ssLouvain
import pomegranate as pmg

from icat.src import utils

class SignalNoiseModel(object):
    def __init__(self):
        """
        Model technical noise in a single-cell expression profile.
        
        Attributes
        ----------
        gmm_ : pomegranate.GeneralMixtureModel 
            Mixture of a Poisson and Normal Distribution to de-convolve noise
            and signal.
        range_ : list
            Two element list with range of possible values. Minimum is always
            zero, and max is set to the maximum observed value + spread.
        weights_ : dict
            Mixing proportions beteen technical noise and actual signal. Keyed
            by 'noise' and 'signal', respectively.

        Methods
        -------
        pdf : Calculate the probability of values within an array.
        cdf : Calculate cumulative density probabilities of values in an array.
        threshold: Find the first value such that P(Noise) < P(Signal)
        """
        self.gmm_ = None
        self.range_ = None
        self.weights_ = None

    def fit(self, X):
        """
        Fit data to a mixture model to distinguish signal from technical noise.

        Fit data to a mixture model where technical noise of scRNAseq profiles
        is modelled using a Poisson distribution, and true expression is
        modelled as a Gaussian Distribution. 
        
        Parameters
        ----------
        X : numpy.array
            Single-cell expression profile across cells. Values are assumed to
            be log-transformed.
        
        Returns
        -------
        None
        """
        #  use count data for mixtures
        counts, bins = np.histogram(X, bins=30)

        # estimate center as maximum non-zero count
        mode_idx = np.where(counts == np.max(counts[1:]))[0][0]
        # estimate 2SD as center - end value --> find values in normal dist.
        normal_spread = bins[-1] - bins[mode_idx]

        # find minimum value heuristically expected to be in signal
        noise_indices = np.where(bins < bins[mode_idx] - normal_spread)[0]
        if len(noise_indices) > 0:
            normal_min = bins[noise_indices[-1]]
        else:
            # no values below expected threshold, set to first non-zero value
            normal_min = bins[1]

        # estimate Normal distribution from likely normal samples
        signal = pmg.NormalDistribution.from_samples(X[X >= normal_min])
        # estimate Poisson from likely non-normal samples
        pois_samples = X[X < normal_min]
        percent_zeros = sum(X == 0) / len(X) 
        pois_lambda = percent_zeros
        if len(pois_samples) != 0:
            pois_lambda = max(np.mean(pois_samples), percent_zeros)
        noise = pmg.PoissonDistribution(pois_lambda)
        # instantiate and fit mixture to data
        gmm = pmg.GeneralMixtureModel([noise, signal])
        gmm.fit(X)
        self.gmm_ = gmm
        self.range_ = [0, np.max(X) + normal_spread]
        # exponentiate log weights, will likely change depending on
        # pomegrante 
        self.weights_ = {'noise': np.exp(gmm.weights[0]),
                         'signal': np.exp(gmm.weights[1])}
        self.params_ = {'noise': {'lambda': gmm.distributions[0].parameters[0]},
                        'signal': {'mean': gmm.distributions[1].parameters[0],
                                   'var': gmm.distributions[1].parameters[1]}}

    def __check_fit(self):
        if self.gmm_ is None:
            raise NotFittedError("This SignalNoiseModel is not fitted.")

    def pdf(self, X):
        """
        Calculate probability density function of the fit mixture model.
        
        Parameters
        ----------
        X : np.array
            Linespace defining the domain of the mixture model.
        
        Returns
        -------
        np.array
            Array of probabilities along the domain.
        """
        self.__check_fit()
        return self.gmm_.probability(X)

    def cdf(self, X):
        """
        Calculate cumulative density function of the fit mixture model.
        
        Parameters
        ----------
        X : np.array
            Linespace defining the domain of the mixture model.
        
        Returns
        -------
        np.array
            Array of cumulative densities for each provided value.
        """
        self.__check_fit()
        values = np.zeros_like(X)
        for i, x in enumerate(X):
            space = np.arange(0, x + 0.01, 0.01)
            values[i] = np.sum(self.pdf(space)*0.01)
        return values

    def threshold(self):
        self.__check_fit()
        space = np.arange(self.range_[0], self.range_[1], 0.01).reshape(-1, 1)
        p_x = self.gmm_.predict_proba(space)
        idx = 0
        while idx < p_x.shape[0] - 1 and p_x[idx][0] > p_x[idx][1]:
            idx += 1
        return space[idx][0]

class SEG():
    def __init__(self):
        """Class to find stabley expressed genes between datasets."""
        self.index_ = None

    def __check_fit(self):
        if self.index_ is None:
            raise NotFittedError("This SEG is not fitted.")

    def __corrected_w(self, X, mu_min, mu_max):
        w = sum(X==0) / X.shape[0]
        return np.sqrt(w * (np.mean(X) - mu_min) / (mu_max - mu_min))

    def fit(self, datasets, ds_col):
        combined = utils.rbind_adata(datasets)
        idxs = [np.where(combined.obs[ds_col] == x)[0]\
                for x in combined.obs[ds_col].unique()]
        X = combined.X
        mus = np.mean(X, axis=0)
        mu_min = np.min(mus)
        mu_max = np.max(mus)
        metric_array = np.zeros((X.shape[1], 4))
        for i in range(X.shape[1]):
            gene_X = X[:, i]
            mixture = SignalNoiseModel()
            mixture.fit(np.log2(gene_X + 1))
            # high percentage of noise mixing -> unstable
            mixing = mixture.weights_['noise']
            # high variance in signal -> unstable
            signal_var = mixture.params_['signal']['var']
            # high proportion of zeros -> unstable 
            w = self.__corrected_w(gene_X, mu_min, mu_max)
            # high H -> unstable
            H = stats.kruskal(*[gene_X[idx] for idx in idxs]).statistic
            metric_array[i, :] = np.array([mixing, signal_var, w, H])
        # get ranks for each gene
        sorted_ = np.argsort(metric_array, axis=0)
        ranks = np.empty_like(sorted_)
        for i in range(metric_array.shape[1]):
            ranks[sorted_[:, i], i] = np.arange(metric_array.shape[0])
        ranks = ranks / (ranks.shape[0] - 1)
        self.index_ = np.mean(ranks, axis=1)

    def get_stable(self, n_genes=100):
        self.__check_fit()
        # want lowest
        return np.argsort(self.index_)[:n_genes]


class icat():
    """
    Model to Identify Clusters Across Treatments. 
    """

    def __init__(self, clustering='louvain', treatment_col='treatment',
                 ncfs_kws=None, cluster_kws=None, cluster_col=None,
                 weight_threshold=None, neighbor_kws=None,
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
        return self._treatment_col

    @treatment_col.setter
    def treatment_col(self, value):
        self._treatment_col = value

    @property
    def ncfs_kws(self):
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
        return self._weight_threshold
    
    @weight_threshold.setter
    def weight_threshold(self, value):
        if value is None:
            value = 0
        elif not isinstance(value, (float, int, np.float, np.integer)):
            raise ValueError("Expected numerical value for `weight_threshold`."
                             "Received: {}".format(value))
        self._weight_threshold = value

    @property
    def neighbor_kws(self):
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
        return self._cluster_col
    
    @cluster_col.setter
    def cluster_col(self, value):
        if not isinstance(value, (str, int)) and value is not None:
            raise ValueError('Expected integer index or string name for '+\
                             '`cluster_col`')
        self._cluster_col = value

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
            controls.obs[self.treatment_col] = 'Controls'
        if not isinstance(perturbed, sc.AnnData):
            if isinstance(perturbed, list):
                if not all([isinstance(x, sc.AnnData) for x in perturbed]):
                    raise ValueError("Expected all perturbed datasets to be "
                                     "sc.AnnData objects.")
                if not all([utils.check_matching_genes(controls, x)\
                for x in perturbed]):
                    raise ValueError("Gene columns do not match between control"
                                     " and perturbed cells.")
                perturbed = utils.rbind_adata(perturbed)
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
        # scale perturbed data using control data
        scaler = preprocessing.MinMaxScaler()
        scaler.fit(controls.X)
        # 
        sc.pp.pca(controls, **self.pca_kws)
        sc.pp.neighbors(controls, **self.neighbor_kws)
        sc.tl.umap(controls, min_dist=0.0)
        # no previous clustering provided, cluster ya self
        if self.cluster_col is None:
            if self.clustering == 'louvain':
                sc.tl.louvain(controls, **self.cluster_kws)
                self.cluster_col = 'louvain'
            elif self.clustering == 'leiden':
                sc.tl.leiden(controls, **self.cluster_kws)
                self.cluster_col = 'leiden'
        else:
            if self.cluster_col not in controls.obs.columns:
                raise ValueError("`cluster_col` - {} not found in control data."\
                                 .format(self.cluster_col))
            if np.any([pd.isnull(x) for x in controls.obs[self.cluster_col]]):
                raise ValueError("Expected labelled cells by passing "\
                                 "`cluster_col`={}. ".format(self.cluster_col) +
                                 "Received atleast one unannotated cell.")
        # instantiate ncfs model
        model = NCFS.NCFS(**self.ncfs_kws)

        # scale genes between 0 and 1 -- TODO: change to option 
        fit_X = scaler.transform(controls.X).astype(np.float64)
        # combine control and perturbed data
        combined = utils.rbind_adata([controls, perturbed])
        # fit gene weights using control dataset
        model.fit(fit_X, np.array(controls.obs[self.cluster_col].values))
        # apply learned weights across gene expression matrix
        combined.X = model.transform(combined.X)
        # save genes weights
        combined.var['ncfs.weights'] = model.coef_
        # subset to most informative genes
        selected = combined.var.index.values[
                       np.where(model.coef_**2 > self.weight_threshold)[0]]
        if len(selected) == 0:
            print('WARNING: No feature weights met threshold criteria. '
                  'All genes will be used. Try lowering threshold value for'
                  ' future runs.')
            selected = combined.var.index.values
        combined = combined[:, selected].copy() # copy to be safe
        # create neighbor graph for control+perturbed combined data -- perform
        # pca if more than n_components remain in selected data
        if self.pca_kws['n_comps'] < combined.shape[1]:
            sc.pp.pca(combined, **self.pca_kws)
        sc.pp.neighbors(combined, **self.neighbor_kws)
        sc.tl.umap(combined, min_dist=0.0)
        # grab connectivities of cells
        A_ = combined.uns['neighbors']['connectivities']
        # instantiate semi-supervised Louvain model
        ss_model = ssLouvain.ssLouvain(**self.sslouvain_kws)
        y_ = combined.obs[self.cluster_col]
        # cluster cells
        ss_model.fit(A_, y_)
        # store new cluster labels in cell metadata
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
    model = icat(ncfs_kws={'reg': 3, 'sigma': 2},
                 weight_threshold=0,
                 neighbor_kws={'n_neighbors': 100},
                 sslouvain_kws={'immutable': True, 'precluster': False})
    out = model.cluster(controls, perturbed)
