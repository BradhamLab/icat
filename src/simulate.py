import numpy as np
from scipy import stats


# Don't include batches in initial dataset creation
# create a "pertubation" method that will generate perturbed
# dataset from original.
class SingleCellDataset():
    
    def __init__(self, samples=200, genes=1000, batches=2, 
                 batch_sizes=None, populations=None, pop_sizes=None):
        self.samples = samples
        self.genes = genes
        self.batches = batches
        self.batch_sizes = batch_sizes
        self.populations = populations
        self.pop_sizes = pop_sizes

    @property
    def samples(self):
        return self.samples

    @samples.setter
    def samples(self, value):
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for sample size.")
        self.samples = value

    @property
    def genes(self):
        return self.genes

    @genes.setter
    def genes(self, value):
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for number of genes.")
        self.genes = value

    @property
    def batches(self):
        return self.batches
    
    @batches.setter
    def batches(self, value):
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for number of batches.")

    @property
    def batch_sizes(self):
        return self.batch_sizes
    
    @batch_sizes.setter
    def batch_sizes(self, value):
        if value is None:
            value = np.array([self.samples // self.batches\
                                    for each in range(self.batches)])
            remainder = self.samples % self.batches
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
        if value.sum() != self.samples:
            raise ValueError("Expected batch sizes to total to number of "
                             "samples. Got {}".format(value.sum()))
        self.batch_sizes = value

def average_exp(beta, scale_factor, n=1):
    return beta.rvs(n) * scale_factor

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
    beta = stats.beta(a=2, b=5)
    mus = average_exp(beta, 100, n_genes)
    exp_matrix = simulate_counts(n_samples, mus, beta_0=beta_0)
    return exp_matrix, mus

def simulate_new_batch(n_samples, reference_mus, percentage, beta_0=-1.5):
    effected_genes = np.random.choice(np.arange(reference_mus.size),
                                      size=int(reference_mus.size*percentage))
    exp_shifts = stats.gamma(a=1, scale=1).rvs(size=effected_genes.size)
    reference_mus[effected_genes] *= exp_shifts
    return simulate_counts(n_samples, reference_mus, beta_0)




