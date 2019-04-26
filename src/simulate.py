import numpy as np
from scipy import stats
import pandas as pd
from scanpy import api as sc


# Don't include batches in initial dataset creation
# create a "pertubation" method that will generate perturbed
# dataset from original.
class SingleCellDataset():
    
    def __init__(self, samples=200, genes=1000, populations=2,
                 pop_sizes=None, p_marker=None, dispersion=4, fc=2, scalar=100):
        self.samples = samples
        self.genes = genes
        self.populations = populations
        self.pop_sizes = pop_sizes
        self.p_marker = p_marker
        self.dispersion = dispersion
        self.fc = fc
        self.scalar = scalar

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, value):
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for sample size.")
        self._samples = value

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, value):
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for number of genes.")
        self._genes = value

    @property
    def populations(self):
        return self._populations
    
    @populations.setter
    def populations(self, value):
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for number of "
                             "populations.")
        self._populations = value

    @property
    def pop_sizes(self):
        return self._pop_sizes
    
    @pop_sizes.setter
    def pop_sizes(self, value):
        if value is None:
            value = np.array([self.samples // self.populations\
                                    for each in range(self.populations)])
            remainder = self.samples % self.populations
            if remainder != 0:
                value[:remainder] += 1
        try:
            value = np.array(value, dtype=int)
        except TypeError:
            raise ValueError('Expected numpy.array castable when setting batch '
                             'sizes.')
        except ValueError:
            raise ValueError('Expected numpy.array castable when setting batch '
                             'sizes.')
        if len(value) != self.populations:
            raise ValueError("Number of populations does not match number of "
                             "populations passed. "
                             "Populations set to {}".format(self.populations),
                             ", but received "
                             "{} population sizes.".format(self.pop_sizes))
        if value.sum() != self.samples:
            raise ValueError("Expected batch sizes to total to number of "
                             "samples. Got {}".format(value.sum()))
        self._pop_sizes = value

    @property
    def p_marker(self):
        return self._p_marker

    @p_marker.setter
    def p_marker(self, value):
        if value is None:
            value = 10 / self.genes
        elif value < 0 or value > 1:
            raise ValueError("Expected value for p_marker between 0 and 1. "
                             "Got {}.".format(value))
        self._p_marker = value

    @property
    def dispersion(self):
        return self._dispersion

    @dispersion.setter
    def dispersion(self, value):
        if not isinstance(value, (int, np.integer)):
            raise ValueError("Expected integer value for dispersion parameter "
                             "Received: {} - {}".format(value, type(value)))
        self._dispersion = value

    @property
    def fc(self):
        return self._fc
    
    @fc.setter
    def fc(self, value):
        if not isinstance(value, (float, int, np.float, np.integer)):
            raise ValueError("Expected numerical value for `fc` parameter."
                             "Received: {}".format(type(value)))
        if self.dispersion / value != self.dispersion // value:
            print("Warning: `dispersion` parameter value not cleanly divisible "
                  "by `fc`. Downregulated genes will have dispersions forced to"
                  " nearest integer value.")
        if self.dispersion * value != int(self.dispersion * value):
            print("Warning: float `fc` parameter will result in float "
                  "dispersion value for upregulated genes: upregulated value "
                  "will be set to nearest integer.")
        self._fc = value

    @property
    def scalar(self):
        return self._scalar

    @scalar.setter
    def scalar(self, value):
        if not isinstance(value, (float, int, np.float, np.integer)):
            raise ValueError("Expected numerical value for `scalar` parameter."
                             "Received: {}".format(type(value)))
        if value < 1:
            raise ValueError("Expected `scalar` value > 1: values less than one"
                             " will result in average gene expressions < 1.")
        self._scalar = value

    def __repr__(self):
        header = "Simulated Single-Cell Dataset\n"
        out = ["{}: {}".format(k, v) for k, v in self.get_params().items()]
        return header + '\n'.join(out)

    def get_params(self):
        out_dict = {'samples': self.samples,
                    'genes': self.genes,
                    'populations': self.populations,
                    'pop_sizes': self.pop_sizes,
                    'p_marker': self.p_marker}
        return out_dict

    def simulate(self):
        X_ = np.zeros((self.samples, self.genes), dtype=int)
        mus_ = average_exp(scale_factor=self.scalar, n=self.genes)
        disp_ = np.ones(mus_.size)*2
        percentiles = np.percentile(mus_, [25, 75])
        down_r = self.dispersion // self.fc
        up_r = int(self.dispersion * self.fc)
        obs = pd.DataFrame(index=range(self.samples), columns=["Population"])
        var = pd.DataFrame(index=range(self.genes),
                           columns=['Pop.{}.Marker'.format(i + 1) for i in\
                                                       range(self.populations)])
        var.fillna(False, inplace=True)
        for i in range(self.populations):
            var['Pop.{}.Dispersion'.format(i + 1)] = self.dispersion
        var['Base.Mu'] = mus_

        for i in range(self.populations):
            n_markers = stats.binom(self.genes, self.p_marker).rvs()
            markers = np.random.choice(np.arange(self.genes), n_markers)
            shifts = np.zeros(shape=n_markers)
            for k, j in enumerate(markers):
                # gene is near medial expression values, up or down regulate
                if percentiles[0] <= mus_[j] <= percentiles[1]:
                    shifts[k] = np.random.choice([down_r, up_r], size=1)
                # gene is below medial expression values, up regulate
                elif mus_[j] < percentiles[0]:
                    shifts[k] = up_r
                # gene is above medial expression values, down regulate
                else:
                    shifts[k] = down_r
            # shift marker gene expression values away from baseline averages.
            pop_disp_ = disp_.copy()
            pop_disp_[markers] = shifts
            if i == 0:
                start = 0
            else:
                start = self.pop_sizes[:i].sum()
            X_[start:start + self.pop_sizes[i], :] = simulate_counts(
                                                              self.pop_sizes[i],
                                                              mus_, r=pop_disp_)
            obs.loc[start:start + self.pop_sizes[i], 'Population'] = i + 1
            var.loc[markers, 'Pop.{}.Marker'.format(i + 1)] = True
            var.loc[:, 'Pop.{}.Dispersion'.format(i + 1)] = pop_disp_
        return sc.AnnData(X=X_, obs=obs, var=var)

# TODO: add option to simulate untargetted populations/ add treatment
# specific populations
def perturb(andata, samples=200, pop_targets=None, gene_targets=None,
            percent_perturb=None, pop_sizes=None):
    # check pop_targets in andata
    if pop_targets is None:
        pop_targets = andata.obs['Population'].unique()
        pop_ratios = andata.obs['Population'].value_counts().values
        pop_ratios = pop_ratios / pop_ratios.sum()
    else:
        pop_targets = __check_np_castable(pop_targets, 'pop_targets')
        if not set(pop_targets).issubset(andata.obs['Population']):
            diff = set(pop_targets).difference(andata.obs['Population'])
            raise ValueError("Undefined populations: {}".format(diff))
    # check population sizes, if none, match with ratio in andata
    if pop_sizes is None:
        try:
            pop_ratios
        except NameError:
            pop_ratios = np.ones(pop_targets.size) * 1 / len(pop_targets)
        pop_sizes = (pop_ratios * samples).astype(int)
        print(pop_sizes)
        if samples % pop_sizes.sum() != 0:
            remainder = samples % pop_sizes.sum()
            iters = 0
            while remainder != 0:
                pop_sizes[iters % len(pop_sizes)] += 1
                remainder = samples % pop_sizes.sum()
    else:
        pop_sizes = __check_np_castable(pop_sizes, 'pop_sizes')
        if pop_sizes.sum() != samples:
            raise ValueError('Population sizes do not sum to number of samples.')
    # determine which genes to perturb
    if gene_targets is not None:
        if not set(gene_targets).issubset(range(andata.shape[1])):
            raise ValueError("Unrecognized gene targets: {}".format(
                          set(gene_targets).difference(range(andata.shape[1]))))
        gene_targets = __check_np_castable(gene_targets, 'gene_targets')
            
    if percent_perturb is not None and gene_targets is not None:
        n_genes = int(percent_perturb * andata.shape[1])
        if len(gene_targets) < n_genes:
            n_genes -= len(gene_targets)
            p_genes = set(range(andata.shape[1])).difference(gene_targets)
            targets = np.random.choice(p_genes, n_genes)
            gene_targets = np.hstack((targets, gene_targets))
    else:
        if percent_perturb is None:
            percent_perturb = 0.2
        gene_targets = np.random.choice(int(andata.shape[1] * percent_perturb),
                                        andata.shape[1])

    X_ = np.zeros((samples, andata.shape[1]))
    mus_ = andata.var['Base.Mu'].values
    exp_shifts = stats.gamma(a=1, scale=1).rvs(size=gene_targets.size)
    populations = []
    for i, each in enumerate(pop_targets):
        populations += [each] * pop_sizes[i]
    obs_ = pd.DataFrame(populations, columns=['Population'])
    var_ = andata.var.copy()
    var_['Perturbation.Shift'] = None
    var_.loc[gene_targets, 'Perturbation.Shift'] = exp_shifts
    for i, each in enumerate(pop_targets):
        pop_disp = andata.var['Pop.{}.Dispersion'.format(each)].values
        pop_mus = mus_ * exp_shifts
        if i == 0:
            start = 0
        else:
            start = pop_sizes[i - 1]
        X_[start:start + pop_sizes[i]] = simulate_counts(pop_sizes[i],
                                                         pop_mus, r=pop_disp)
    return sc.AnnData(X=X_, obs=obs_, var=var_)


def __check_np_castable(obj, name):
    if not isinstance(obj, np.ndarray):
        try:
            obj = np.array(obj)
        except:
            raise ValueError("Expected numpy.ndarray castable object for "
                            "`{}`. Got {}.".format(obj, type(obj)))
    return obj


def average_exp(scale_factor, n=1):
    return stats.beta(a=2, b=5).rvs(n) * scale_factor


def sigmoid(x):
    return 1 / (1 + np.e ** -x)


def dropout_probability(mu, median_avg, beta_0=-1.5):
    x = beta_0 + 1 / median_avg * mu
    return sigmoid(x)


def sample_count(mu, p, n=200, r=2):
    dropout = stats.bernoulli(p).rvs(n)
    counts = stats.nbinom(r, 1 - mu / (mu + 1)).rvs(n)
    # if dropout == 0, gene dropped out, multiplying goes to zero
    return counts * dropout


def simulate_counts(n_samples, mus, r=2, beta_0=-1.5):
    exp_matrix = np.zeros((n_samples, mus.size))
    means = r * mus
    median = np.median(means)
    p_dropout = dropout_probability(means, median, beta_0=beta_0)
    if isinstance(r, int):
        for i in range(means.size):
            exp_matrix[:, i] = sample_count(means[i], p_dropout[i], n_samples, r)
    elif isinstance(r, np.ndarray) and len(r) == len(means):
        for i in range(means.size):
            exp_matrix[:, i] = sample_count(means[i], p_dropout[i], n_samples,
                                            r[i])
    return exp_matrix


def simulate_reference_batch(n_samples=200, n_genes=1000, beta_0=-1.5):
    mus = average_exp(100, n_genes)
    exp_matrix = simulate_counts(n_samples, mus, beta_0=beta_0)
    return exp_matrix, mus


def shift_expression(expression_avgs):
    shifts = stats.gamma(a=1, scale=1).rvs(size=expression_avgs.size)
    return expression_avgs * shifts


def simulate_new_batch(n_samples, reference_mus, percentage, beta_0=-1.5):
    effected_genes = np.random.choice(np.arange(reference_mus.size),
                                      size=int(reference_mus.size*percentage))
    exp_shifts = stats.gamma(a=1, scale=1).rvs(size=effected_genes.size)
    reference_mus[effected_genes] *= exp_shifts
    return simulate_counts(n_samples, reference_mus, beta_0)