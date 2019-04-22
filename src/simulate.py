import numpy as np
from scipy import stats


# Don't include batches in initial dataset creation
# create a "pertubation" method that will generate perturbed
# dataset from original.
class SingleCellDataset():
    
    def __init__(self, samples=200, genes=1000, populations=1,
                 pop_sizes=None, p_marker=None):
        self.samples = samples
        self.genes = genes
        self.populations = populations
        self.pop_sizes = pop_sizes
        self.p_marker = p_marker

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
            value = 7 / self.genes
        elif value < 0 or value > 1:
            raise ValueError("Expected value for p_marker between 0 and 1. "
                             "Got {}.".format(value))
        self._p_marker = value
    
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
        self.X_ = np.zeros((self.samples, self.genes), dtype=int)
        self.markers_ = [[]] * self.populations
        self.marker_mus_ = self.markers_.copy()
        self.mus_ = average_exp(100, n=self.genes)
        for i in range(self.populations):
            n_markers = stats.binom(self.genes, self.p_marker).rvs()
            self.markers_[i] = np.random.choice(np.arange(self.genes),
                                                n_markers)
            self.marker_mus_[i] = shift_expression(self.mus_[self.markers_[i]])
            pop_avgs = self.mus_.copy()
            pop_avgs[self.markers_[i]] = self.marker_mus_[i]
            if i == 0:
                start = 0
            else:
                start = self.pop_sizes[:i].sum()
            print(start, self.pop_sizes[i])
            self.X_[start:start + self.pop_sizes[i], :] = \
                simulate_counts(self.pop_sizes[i], pop_avgs)
            

def average_exp(scale_factor, n=1):
    return stats.beta(a=2, b=5).rvs(n) * scale_factor

def sigmoid(x):
    return 1 / (1 + np.e ** -x)

def dropout_probability(mu, median_avg, beta_0=-1.5):
    x = beta_0 + 1 / median_avg * mu
    return sigmoid(-x)

def sample_count(mu, p, n=200):
    dropout = stats.bernoulli(p).rvs(n)
    counts = stats.nbinom(1, 1 - mu / (mu + 1)).rvs(n)
    # if dropout == 0, gene dropped out, multiplying goes to zero
    return counts * dropout

def simulate_counts(n_samples, mus, beta_0=-1.5):
    exp_matrix = np.zeros((n_samples, mus.size))
    median = np.median(mus)
    p_dropout = dropout_probability(mus, median, beta_0=beta_0)
    for i in range(mus.size):
        exp_matrix[:, i] = sample_count(mus[i], p_dropout[i], n_samples)
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




