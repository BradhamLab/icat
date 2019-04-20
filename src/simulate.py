import numpy as np
from scipy import stats


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

def simulate_reference_batch(n_samples=200, n_genes=1000, beta_0=-1.5):
    exp_matrix = np.zeros((n_samples, n_genes))
    beta = stats.beta(a=2, b=5)
    mus = average_exp(beta, 100, n_genes)
    median = np.median(mus)
    p_dropout = dropout_probability(mus, median, beta_0=beta_0)
    for i in range(n_genes):
        exp_matrix[:, i] = sample_count(mus[i], p_dropout[i], n_samples)
    return exp_matrix, mus



