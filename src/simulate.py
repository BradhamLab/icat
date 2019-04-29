"""
Module to simulate single-cell RNA sequencing data

@author: Dakota Y. Hawkins
@contact: dyh0110@bu.edu 
"""
import numpy as np
from scipy import stats
import pandas as pd
from scanpy import api as sc

class SingleCellDataset():
    """
    Class to simulate a single-cell RNA sequencing dataset.
    
    Attributes
    ----------
    samples : int
        Number of cells to simulate.
    genes : int
        Number of genes to simulate.
    populations : int
        Number of populations to simulate in the dataset.
    pop_sizes : list-like
        Number of cells in each population.
    p_marker : float
        Probability of selecting a gene as a marker gene for a given population.
    dispersion : int, list-like
        Dispersion parameter used to simulate counts in a negative binomial
        model. Either a single, integer value shared across all genes, or an
        array of gene-specific values.
    scalar : float
        Scalar value multiplied to a beta-distributed random variable to
        estimate average gene expression for a simualted gene.
    Methods
    -------
    
    get_params():
        Get model parameters used to simulate dataset.
    simulate():
        Simulate a single-cell RNA sequencing dataset with the set parameters.
    """
    def __init__(self, samples=200, genes=1000, populations=2,
                 pop_sizes=None, p_marker=None, dispersion=1, scalar=100):
        """
        Parameters
        ----------
        samples : int, optional
            Number of cells to simulate. Default is 200.
        genes : int, optional
            Number of genes to simulate. Default is 1000.
        populations : int, optional
            Number of populations to simulate in the dataset. Default is 2.
        pop_sizes : list-like, optional
            Number of cells in each population. Default is None, and the number
            of cells will be split uniformly between populations.
        p_marker : float, optional
            Probability of selecting a gene as a marker gene for a given
            population. Default is None, and will be set so that each
            population will have an average of 10 marker genes.
        dispersion : int, list-like, optional
            Dispersion parameter used to simulate counts in a negative binomial
            model. Either a single, integer value shared across all genes, or an
            array of gene-specific values. Default is 1, and will be shared
            across all genes.
        scalar : float, optional
            Scalar value multiplied to a beta-distributed random variable to
            estimate average gene expression for a simualted gene. Default is
            100.
        """
        self.samples = samples
        self.genes = genes
        self.populations = populations
        self.pop_sizes = pop_sizes
        self.p_marker = p_marker
        self.dispersion = dispersion
        self.scalar = scalar
        
    @property
    def samples(self):
        """Get samples attribute."""
        return self._samples

    @samples.setter
    def samples(self, value):
        """Set samples attribute."""
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for sample size.")
        self._samples = value

    @property
    def genes(self):
        """Get genes attribute."""
        return self._genes

    @genes.setter
    def genes(self, value):
        """Set genes attribute."""
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for number of genes.")
        self._genes = value

    @property
    def populations(self):
        """Get populations attribute."""
        return self._populations
    
    @populations.setter
    def populations(self, value):
        """Set populations attribute."""
        if value < 0 and int(value) != value:
            raise ValueError("Expected positive integer for number of "
                             "populations.")
        self._populations = value

    @property
    def pop_sizes(self):
        """Get pop_sizes attribute."""
        return self._pop_sizes
    
    @pop_sizes.setter
    def pop_sizes(self, value):
        """Set pop_sizes attribute."""
        if value is None:
            value = np.array([self.samples // self.populations\
                                    for each in range(self.populations)])
            remainder = self.samples % self.populations
            if remainder != 0:
                value[:remainder] += 1
        try:
            value = np.array(value, dtype=int)
        except TypeError:
            raise ValueError('Expected numpy.array castable when setting '
                             'population sizes.')
        except ValueError:
            raise ValueError('Expected numpy.array castable when setting '
                             'population sizes.')
        if len(value) != self.populations:
            raise ValueError("Number of populations does not match number of "
                             "populations passed. "
                             "Populations set to {}".format(self.populations),
                             ", but received "
                             "{} population sizes.".format(self.pop_sizes))
        if value.sum() != self.samples:
            raise ValueError("Expected population sizes to total to number of "
                             "samples. Got {}".format(value.sum()))
        self._pop_sizes = value

    @property
    def p_marker(self):
        """Get p_marker attribute."""
        return self._p_marker

    @p_marker.setter
    def p_marker(self, value):
        """Set p_marker attribute."""
        if value is None:
            value = 10 / self.genes
        elif value < 0 or value > 1:
            raise ValueError("Expected value for p_marker between 0 and 1. "
                             "Got {}.".format(value))
        self._p_marker = value

    @property
    def dispersion(self):
        """Get dispersion parameter."""
        return self._dispersion

    @dispersion.setter
    def dispersion(self, value):
        """Set dispersion parameter."""
        if isinstance(value, (list, np.ndarray)):
            if not len(value) == self.genes:
                raise ValueError("Number of dispersions passed does not match "
                                 "the number of genes passed. Dispersion should"
                                 " either be a list-like of length `genes` or "
                                 "a single integer value.")
            value = np.array(value, dtype=int)
            if not np.all(value > 0):
                raise ValueError("Expected non-negative positive integer for "
                                 "all dispersion parameters.")
        elif not isinstance(value, (int, np.integer)):
            raise ValueError("Expected integer value for dispersion parameter "
                             "Received: {} - {}".format(value, type(value)))
            if value < 1:
                raise ValueError("Dispersion must be a non-negative integer. "
                                 "Received: {}".format(value))
        self._dispersion = value

    @property
    def scalar(self):
        """Get scalar parameter."""
        return self._scalar

    @scalar.setter
    def scalar(self, value):
        """Set scalar parameter."""
        if not isinstance(value, (float, int, np.float, np.integer)):
            raise ValueError("Expected numerical value for `scalar` parameter."
                             "Received: {}".format(type(value)))
        if value < 1:
            raise ValueError("Expected `scalar` value > 1: values less than one"
                             " will result in average gene expressions < 1.")
        self._scalar = value

    def __repr__(self):
        """Return string representation of SingleCellDataset object."""
        header = "Simulated Single-Cell Dataset\n"
        out = ["{}: {}".format(k, v) for k, v in self.get_params().items()]
        return header + '\n'.join(out)

    def get_params(self):
        """Get parameter set used to simulate dataset."""
        out_dict = {'samples': self.samples,
                    'genes': self.genes,
                    'populations': self.populations,
                    'pop_sizes': self.pop_sizes,
                    'p_marker': self.p_marker,
                    'dispersion': self.dispersion,
                    'scalar': self.scalar}
        return out_dict

    def simulate(self):
        """
        Simulate read counts across genes and populations.           
        
        Parameters
        ----------
        None
        
        Returns
        -------
        sc.AnnData
            Annotated dataframe of simulated data.

        Notes
        -----
        Counts are simulated from a negative binomial distribution with dropout.

        ..math::
             c_{ij} ~ NV(r_ij, 1 - \dfrac{\mu_{ij}}{\mu + 1} | Bernoulli(p_ij))

        Where :math:`c_ij` is the read counts for gene :math:`j` in cell
        :math:`i`, :math:`r_ij` is the dispersion parameter for gene `j` when
        expressed in cell `i`, :math:`\mu_{ij}` is the average expression value
        for gene :math:`j` in cell :math:`i`, and :math:`p_ij` is the dropout
        probability of gene :math:`j` in cell :math:`i`. As :math:`c_{ij}` is
        conditioned on :math:`p_{ij}`, if :math:`\not p_{ij}, c_{ij} = 0`. 

        The baseline average expression value for gene `j`, :math:`\mu_{j}`, is
        modelled using a scaled Beta distribution with parameters :math:`a = 2`
        and :math:`b = 5`. Each Beta random variable is scaled by a constant
        :math:`c = 100` to calculate the final average.
        
        ..math::
            \mu_j \sim Beta(a, b) \cdot c

        If gene :math:`j` was selected as a marker gene for population
        :math:`K`, this average is scaled again to represent up or down
        regulation in population :math:`K` compared to baseline. The marker gene
        scalar is drawn from a Gamma-distributed random variable, where:

        ..math::
            \gamma \sim Gamma(\alpha, \beta) \\
            \mu_{j}^K = \gamma \mu_j \\
            
        and :math:`\alpha = 2`, `\beta = 2`. Thus, letting :math:`M` represent
        the set of all marker genes.
            
        ..math::
            \mu_ij = 
            \begin{cases}
                \mu_{j}^K & \forall i \in K \and j \in M`, \\
                \mu_j & otherwise
            \end{cases}

        The probability for a dropout event for gene :math:`j` in cell `i` is
        modelled by passing an affine transformation through a sigmoid function:

        ..math::
            p_ij = sigmoid(x) \\
            sigmoid(x) = \dfrac{1}{1 + e^{-x}} \\
            x_ij = \beta_0 + \dfrac{mu_ij}{median(\vec \mu)} \\
            \vec \mu = \{\mu_11, \mu_12, \ldots , \mu_n{p - 1}, \mu_np}
        
        Here, :math:`beta_0 = -1.5` and :math:`median(\vec \mu)` is the median
        average expression value across all genes over all cells.
        """
        X_ = np.zeros((self.samples, self.genes), dtype=int)
        mus_ = average_exp(scale_factor=self.scalar, n=self.genes)
        obs = pd.DataFrame(index=range(self.samples), columns=["Population"])
        var = pd.DataFrame(index=range(self.genes),
                           columns=['Pop.{}.Marker'.format(i + 1) for i in\
                                                       range(self.populations)])
        var.fillna(False, inplace=True)
        for i in range(self.populations):
            var['Pop.{}.Mu'.format(i + 1)] = mus_
        var['Base.Dispersion'] = self.dispersion
        var['Base.Mu'] = mus_
        gamma = stats.gamma(a=2, scale=2)
        for i in range(self.populations):
            # randomly select number of marker genes, must have at least 1
            n_markers = max(stats.binom(self.genes, self.p_marker).rvs(), 1) 
            markers = np.random.choice(np.arange(self.genes), n_markers)
            shifts = gamma.rvs(n_markers)
            # shift marker gene expression values away from baseline averages.
            pop_mus_ = mus_.copy()
            pop_mus_[markers] *= shifts
            if i == 0:
                start = 0
            else:
                start = self.pop_sizes[:i].sum()
            X_[start:start + self.pop_sizes[i], :] = simulate_counts(
                                                              self.pop_sizes[i],
                                                              pop_mus_,
                                                              r=self.dispersion)
            obs.loc[start:start + self.pop_sizes[i], 'Population'] = i + 1
            var.loc[markers, 'Pop.{}.Marker'.format(i + 1)] = True
            var.loc[:, 'Pop.{}.Mu'.format(i + 1)] = pop_mus_
        return sc.AnnData(X=X_, obs=obs, var=var)

# TODO: add option to simulate untargetted populations/ add treatment
# specific populations
def perturb(andata, samples=200, pop_targets=None, gene_targets=None,
            percent_perturb=None, pop_sizes=None):
    """
    Perturb a simulated single-cell dataset.

    Perturbs a simulated single-cell dataset by applying a Gamma-distributed
    scalar shift to the average expression value to a set number of genes.

    That is, for a specific simulated gene following a negative binomial 
    distribution :math:`gene_i \sim NB(r_i, \mu_i)`, the perturbation is
    modelled using the following equations: 

    .. math::
        \gamma_i \sim Gamma(2, 2) \\
        gene_{i} \sim NB(r_i, \gamma_i \cdot \mu_i)

    Where :math:`mu_i` is taken from the reference dataset.

    Parameters
    ----------
    andata : sc.AnnData
        Simulated annotated dataframe to be perturbed. Output from
        SingleCellDataset.simulate().
    samples : int, optional
        Number of perturbed cells to simulate, by default 200.
    pop_targets : list-like, optional
        Populations to simulate, by default None, and all populations present in
        `andata` will be simulated.
    gene_targets : list-like, optional
        Genes to perturb, by default None, and targets will be randomly chosen.
    percent_perturb : float, optional
        Percentage of genes to perturb. By default None, and if no argument is
        provided for `gene_targets`, will be set to 20%. If *both* arguments
        for `gene_targets` and `percent_perturb` are provided, the remainder of
        genes between `gene_targets` and the number of genes dictated by
        `percent_purturbed`, will be randomly selected.
    pop_sizes : list-like, optional
        Desired population sizes for each perturbed population. By default None,
        and the number of cells will be evenly distributed between populations.
    
    Returns
    -------
    sc.AnnData
        Perturbed single-cell dataset.
    
    Raises
    ------
    ValueError
        Raised if a list-like argument is not castable to a numpy.ndarray.
    ValueError
        Raised target populations are not contained in provied annotated
        dataframe
    ValueError
        Raised if provided population sizes do not sum to number of samples
        passed.
    ValueError
        Raised if provided gene targets are not contained in the provided
        annotated dataframe.
    """
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
    disp_ = andata.var['Base.Dispersion'].values
    exp_shifts = np.ones(andata.shape[1])
    exp_shifts[gene_targets] = stats.gamma(a=2, scale=2).\
                                     rvs(size=gene_targets.size)
    populations = []
    for i, each in enumerate(pop_targets):
        populations += [each] * pop_sizes[i]
    obs_ = pd.DataFrame(populations, columns=['Population'])
    var_ = andata.var.copy()
    var_['Perturbation.Shift'] = exp_shifts
    for i, each in enumerate(pop_targets):
        pop_mus = andata.var['Pop.{}.Mu'.format(each)].values * exp_shifts
        if i == 0:
            start = 0
        else:
            start = pop_sizes[:i].sum()
        X_[start:start + pop_sizes[i]] = simulate_counts(pop_sizes[i],
                                                         pop_mus, r=disp_)
    return sc.AnnData(X=X_, obs=obs_, var=var_)


def __check_np_castable(obj, name):
    """Check whether an object is castable to a numpy.ndarray."""
    if not isinstance(obj, np.ndarray):
        try:
            obj = np.array(obj)
        except:
            raise ValueError("Expected numpy.ndarray castable object for "
                            "`{}`. Got {}.".format(obj, type(obj)))
    return obj


def average_exp(scale_factor, n=1):
    """
    Simulate average expression parameters for simulated genes.

    Simulate average expression parameters for simulated genes. Average values
    are modelled using a scaled beta distribution.

    ..math:
        \mu_i \sim Beta(a, b) \cdot c

    Where :math:`a = 2`, :math:`b = 5`, and :math:`c` is a user provided scale
    factor. 
    
    Parameters
    ----------
    scale_factor : float
        Scalar to multiply beta-distributed random variable to calculate
        expression averages. 
    n : int, optional
        Number of genes to model, by default 1. 
    
    Returns
    -------
    np.ndarray
        Average gene expression values. 
    """
    return stats.beta(a=2, b=5).rvs(n) * scale_factor


def sigmoid(x):
    """
    Sigmoid function used to estimate probability of dropout.
    
    Parameters
    ----------
    x : float
        Value to pass to sigmoid function.
    
    Returns
    -------
    float
        Value of `x` passed through sigmoid function.        
    """
    return 1 / (1 + np.e ** -x)


def dropout_probability(mu, median_avg, beta_0=-1.5):
    """
    Estimate the probability of dropout for a given gene.

    Estimate the probability of a dropout even using a sigmoid function, and
    the following models. 

    ..math::
        p_i = sigmoid(x) \\
        sigmoid(x) = \dfrac{1}{1 + e^{-x}} \\
        x_i = \beta_0 + \dfrac{mu_i}{median(\vec \mu)} \\
        \vec \mu = \{\mu_1, \mu_2, \ldots , \mu_{p - 1}, \mu_p}
    
    Parameters
    ----------
    mu : float, numpy.ndarray
        Either the average expression value for a single gene, or an array of
        average expression values for a set of genes.
    median_avg : float
        Median average expression value over all genes.
    beta_0 : float, optional
        Affine term to add to scalar multiple of `mu`. By default -1.5.
    
    Returns
    -------
    float
        Probability of gene(s) to experience a dropout event.
    """
    x = beta_0 + 1 / median_avg * mu
    return sigmoid(x)

# TODO: change mu + 1 to mu + r
def sample_count(mu, p, n=200, r=2):
    """
    Sample gene counts from a negative binomial distribution with dropout. 

    Sample gene counts from a negative binomial distribution with dropout. Let
    :math:`c_{ij}` represent the counts for cell :math:`i` over gene :math:`j`.
    Then,

    ..math::
        c_{ij} ~ NV(r_ij, 1 - \dfrac{\mu_ij}{mu + 1} | Bernoulli(p_ij))

    If :math:`\not Bernoulli(p_ij)`, `c_ij = 0`. 

    Parameters
    ----------
    mu : float
        Gene expression average for specific cell-gene combo.
    p : float
        Probability of expression dropout for specific cell-gene combo.
    n : int, optional
        Number of cells to simulate, by default 200.
    r : int, optional
        Dispersion paramter for negative binomial model, by default 2.
    
    Returns
    -------
    np.ndarray
        Array of gene counts for each cell.
    """
    dropout = stats.bernoulli(p).rvs(n)
    counts = stats.nbinom(r, 1 - mu / (mu + 1)).rvs(n)
    # if dropout == 0, gene dropped out, multiplying goes to zero
    return counts * dropout


def simulate_counts(n_samples, mus, r=2, beta_0=-1.5):
    """
    Simulate counts across genes for a set number of samples.
    
    Parameters
    ----------
    n_samples : int
        Number of samples to simulate.
    mus : np.ndarray
        Expression averages for each gene. 
    r : int, numpy.ndarray, optional
        Dispersion parameter in negative binomial model. Can either be a single
        integer value shared between all genes, or a numpy array containing
        gene-specific values. By default 2.
    beta_0 : float, optional
        Affine constant used in estimating dropout, by default -1.5.
    
    Returns
    -------
    numpy.ndarray
        An :math:`n \times p` count matrix, where :math:`n` is number of
        samples, defined by `n_samples`, and :math:`p` is the number of genes,
        defined by the size of `mus`. 
    
    Raises
    ------
    ValueError
        Raised if a numpy array is passed for `r` and the lengths of `r` and
        `mus` do not align.
    ValueError
        Raised if `r` is neither a numpy array or an integer value.
    """
    exp_matrix = np.zeros((n_samples, mus.size))
    means = r * mus
    median = np.median(means)
    p_dropout = dropout_probability(means, median, beta_0=beta_0)
    if isinstance(r, (int, np.integer)):
        for i in range(means.size):
            exp_matrix[:, i] = sample_count(means[i], p_dropout[i], n_samples,
                                            r)
    elif isinstance(r, np.ndarray) and len(r) == len(means):
        for i in range(means.size):
            exp_matrix[:, i] = sample_count(means[i], p_dropout[i], n_samples,
                                            r[i])
    else:
        if isinstance(r, np.ndarray):
            raise ValueError("Length of `r` and `mu` do not match.")
        else:
            raise ValueError("Expected integer or numpy array for `r`. "
                             "Received: {}.".format(type(r)))
    return exp_matrix
