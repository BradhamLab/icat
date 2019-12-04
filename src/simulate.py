"""
Module to simulate single-cell RNA sequencing data

@author: Dakota Y. Hawkins
@contact: dyh0110@bu.edu 
"""
import re

import numpy as np
import pandas as pd
from scipy import stats

from icat.src import utils
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
                 pop_sizes=None, p_marker=None, dispersion=1, scalar=100,
                 percentile=50, method='conditional', dropout=0.66):
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
        percentile : float, optional
            Float value between 0 and 100 denoting which percentile to use when
            calculating dropout probabilities. Default is 50, and the median
            will be calculated. 
        """
        self.samples = samples
        self.genes = genes
        self.populations = populations
        self.pop_sizes = pop_sizes
        self.p_marker = p_marker
        self.dispersion = dispersion
        self.scalar = scalar
        self.percentile = percentile
        self.method = method
        self.dropout = dropout
        
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
        elif isinstance(value, (int, np.integer)):
            if value < 1:
                raise ValueError("Dispersion must be a non-negative integer. "
                                 "Received: {}".format(value))
            # cast to array of length `genes`
            value = np.array([value]*self.genes) 
        else:
            raise ValueError("Expected integer value for dispersion parameter "
                             "Received: {} - {}".format(value, type(value)))

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
    
    @property
    def percentile(self):
        """Get percentile parameter."""
        return self._percentile
    
    @percentile.setter
    def percentile(self, value):
        try:
            float(value)
        except:
            raise ValueError("Expected numerical value for `percentile`.")
        if not 0 <= value <= 100:
            raise ValueError("`percentile` must be a float between 0 and 100.")
        self._percentile = value

    @property
    def method(self):
        """
        Which method to use to calculate gene dropout. Options are 
        "conditional" and "uniform".
        """
        return self._method

    @method.setter
    def method(self, value):
        if value not in ['conditional', 'uniform']:
            raise ValueError("Unsupported method {}.".format(value))
        self._method = value

    @property
    def dropout(self):
        """
        Average percent of genes to dropout during count simulation. Only used
        if `method == uniform`.
        """
        return self._dropout
    
    @dropout.setter
    def dropout(self, value):
        if value < 0 or value > 1:
            raise ValueError("Expected values between zero and 1. "
                             "Received {}.".format(value))
        self._dropout = value
        

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
                    'scalar': self.scalar,
                    'percentile': self.percentile,
                    'method': self.method,
                    'dropout': self.dropout}
        return out_dict

    def simulate(self):
        r"""
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
            \vec \mu = \{\mu_1, \mu_1, \ldots , \mu_{p - 1}, \mu_p}
        
        Here, :math:`beta_0 = -1.5` and :math:`median(\vec \mu)` is the median
        average baseline expression value across all genes. That is, shifts
        for marker genes are not taken into consideration.
        """
        # data frame to track cell annotations
        obs = pd.DataFrame(index=range(self.samples),
                           columns=["Population"])
        # data frame to track gene annotations
        var = pd.DataFrame(index=range(self.genes),
                           columns=['Pop.{}.Marker'.format(i + 1) for i in\
                                                       range(self.populations)])
        var.fillna(False, inplace=True)
        var['Base.Dispersion'] = self.dispersion
        # get baseline expression averages
        mus_ = average_exp(scale_factor=self.scalar, n=self.genes)
        var['Base.Mu'] = mus_
        # create gene x pop matrix to hold pop-specific expression averages
        mus = np.ones((mus_.size, self.populations)) * mus_.reshape(-1, 1)
        # get number of marker genes for each population
        n_markers = stats.binom(self.genes, self.p_marker).rvs(self.populations)
        # must have at least 1 marker gene
        n_markers[n_markers == 0] = 1
        # distribution to sample expression shifts from 
        gamma = stats.gamma(a=2, scale=2)
        possible_markers = var.index.values
        for i, n in enumerate(n_markers):
            # pick marker genes and shift their expression in population i
            markers = np.random.choice(possible_markers, n)
            mus[markers, i] = mus[markers, i] * gamma.rvs(n)
            # log marker genes in var data frame
            var.loc[markers, 'Pop.{}.Marker'.format(i + 1)] = True
            # remove selected markers from pool, so that populations don't
            # share marker genes
            possible_markers = np.array(list(
                                    set(possible_markers).difference(markers)))
        # log population averages
        for i in range(self.populations):
            var['Pop.{}.Mu'.format(i + 1)] = mus[:, i]
        
        X, labels = simulate_counts(self.samples, mus, self.dispersion,
                                    self.populations, self.pop_sizes,
                                    percentile=self.percentile,
                                    method=self.method, dropout=self.dropout)
        obs['Population'] = labels
        # rename obs and var rows for clarity
        obs.rename(index={i:'cell-{}'.format(i + 1) for i in obs.index}, inplace=True)
        var.rename(index={i:'gene-{}'.format(i + 1) for i in var.index}, inplace=True)
        return sc.AnnData(X=X, obs=obs, var=var)

def dispersions(size, a=1, b=5):
    """
    Randomly choose dispersion parameter 'r' for simulating cell counts.
    
    Parameters
    ----------
    size : int, tuple
        Size of array to create.
    a : int, optional
        Minimum dispersion value, by default 1.
    b : int, optional
        Maximum dispersion value, by default 5.
    
    Returns
    -------
    numpy.ndarray
        array of dispersion paramters.
    """
    return np.random.choice(range(a, b), size=size)

def population_markers(adata):
    marker_cols = [x for x in adata.var.columns if 'Marker' in x]
    markers = {i + 1: adata.var[adata.var[x]].index.values\
               for i, x in enumerate(marker_cols)}
    return markers

# TODO: add option to simulate untargetted populations/ add treatment
# specific populations
def perturb(adata, samples=200, pop_targets=None, gene_targets=None,
            percent_perturb=None, pop_sizes=None, percentile=50):
    r"""
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
    adata : sc.AnnData
        Simulated annotated dataframe to be perturbed. Output from
        SingleCellDataset.simulate().
    samples : int, optional
        Number of perturbed cells to simulate, by default 200.
    pop_targets : list-like, optional
        Populations to simulate, by default None, and all populations present in
        `adata` will be simulated.
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
    percentile : float, optional
        Percentile to use when calculationg dropout probabilities. Default is
        0.5, and the median will be used.
    
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
    # check pop_targets in adata
    if pop_targets is None:
        pop_targets = adata.obs['Population'].unique()
        pop_ratios = adata.obs['Population'].value_counts().values
        pop_ratios = pop_ratios / pop_ratios.sum()
    else:
        pop_targets = utils.check_np_castable(pop_targets, 'pop_targets')
        if not set(pop_targets).issubset(adata.obs['Population']):
            diff = set(pop_targets).difference(adata.obs['Population'])
            raise ValueError("Undefined populations: {}".format(diff))
    # check population sizes, if none, match with ratio in adata
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
        pop_sizes = utils.check_np_castable(pop_sizes, 'pop_sizes')
        if pop_sizes.sum() != samples:
            raise ValueError('Population sizes do not sum to number of samples.')
    # determine which genes to perturb
    if gene_targets is not None:
        if not set(gene_targets).issubset(adata.var.index):
            raise ValueError("Unrecognized gene targets: {}".format(
                          set(gene_targets).difference(adata.var.index)))
        gene_targets = utils.check_np_castable(gene_targets, 'gene_targets')
            
    if percent_perturb is not None and gene_targets is not None:
        n_genes = int(percent_perturb * adata.shape[1])
        if len(gene_targets) < n_genes:
            n_genes -= len(gene_targets)
            markers = population_markers(adata)
            ignore_genes = gene_targets
            for x in markers.values():
                ignore_genes = np.hstack((ignore_genes, x))
            p_genes = list(set(adata.var.index).difference(ignore_genes))
            targets = np.random.choice(p_genes, n_genes)
            gene_targets = np.hstack((targets, gene_targets))
    elif gene_targets is None:
        if percent_perturb is None:
            percent_perturb = 0.2
        gene_targets = np.random.choice(adata.var.index,
                                        int(adata.shape[1] * percent_perturb))
    markers = population_markers(adata)
    disp_ = adata.var['Base.Dispersion'].values.reshape((-1, 1))
    var_ = adata.var.copy()
    var_['Perturbation.Shift'] = np.ones((adata.shape[1], 1))
    var_.loc[gene_targets, 'Perturbation.Shift'] = stats.gamma(a=2, scale=2).\
                                                   rvs(size=gene_targets.size)
    # log population identity in perturbed data
    # if marker gene is modified, populations is considered perturbed
    populations = []
    for i, each in enumerate(pop_targets):
        markers_i = markers[each]
        if len(set(markers_i).intersection(gene_targets)) != 0:
            name = 'Perturbed.{}'.format(each)
        else:
            name = str(each)
        populations += [name] * pop_sizes[i]
    obs_ = pd.DataFrame(populations, columns=['Population'],
                        index=["cell-{}".format(i + 1) for i in\
                               range(len(populations))])
    pop_columns = ['Pop.{}.Mu'.format(x) for x in pop_targets]
    # calculate control median of averages to ensure equal dropout rates between
    # datasets.
    # control_median_ = np.median(adata.var[pop_columns].values\
    #                 * np.ones((adata.shape[1], len(pop_targets)))\
    #                 * disp_)
    # calculate average expression values for each gene in eac population
    # in perturbed dataset
    mus = adata.var[pop_columns].values \
        * np.ones((adata.shape[1], len(pop_columns))) \
        * var_['Perturbation.Shift'].values.reshape(-1, 1)
    X_, __ = simulate_counts(samples, mus, disp_, len(pop_targets), pop_sizes,
                             percentile=percentile)
    return sc.AnnData(X=X_, obs=obs_, var=var_)


class Experiment(object):
    """Class to simulate scRNA experiments with perturbations."""

    def __init__(self, control_kwargs=None, perturb_kwargs=None):
        """
        Simulate a scRNAseq experiment with perturbations. 
        
        Parameters
        ----------
        control_kwargs : dict, optional
            Dictionary of keyword arguments that specify simulation parameters.
            See `SingleCellDataSet` for more infomration. By default None, and
            default parameters for SingleCellDataSet will be used. 
        perturb_kwargs : dict, optional
            [description], by default None
        """

        self.control_kwargs = control_kwargs
        self.perturb_kwargs = perturb_kwargs

    @property
    def control_kwargs(self):
        """Get control_kwargs attribute."""
        return self._control_kwargs

    @control_kwargs.setter
    def control_kwargs(self, value):
        """Set control_kwargs attribute."""
        default_kws = utils.get_default_kwargs(SingleCellDataset, ['self'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'control_kwargs')
        else:
            value = default_kws
        self._pca_kws = value
        self._control_kwargs = value

    @property
    def perturb_kwargs(self):
        """Get perturb_kwargs attribute."""
        return self._perturb_kwargs

    @perturb_kwargs.setter
    def perturb_kwargs(self, value):
        """Set perturb_kwargs attribute."""
        default_kws = utils.get_default_kwargs(perturb, ['adata'])
        if value is not None:
            value = utils.check_kws(default_kws, value, 'perturb_kwargs')
        else:
            value = default_kws
        self._pca_kws = value
        if 'percentile' not in value and 'percentile' in self.control_kwargs:
            value['percentile'] = self.control_kwargs['percentile']
        self._perturb_kwargs = value

    def simulate_controls(self):
        adata = SingleCellDataset(**self.control_kwargs).simulate()
        return adata

    def new_population(self, adata, n_cells, perturbed=False, pop_id=None):
        """
        Add a new populatin to a previously simulated single-cell dataset.
        
        Parameters
        ----------
        adata : sc.AnnData
            Previously simulated single-cell dataset.
        n_cells : int
            Number of cells to simulate.
        perturbed : bool, optional
            Whether the previously simulated dataset is a perturbed dataset.
            Default is False.
        pop_id : int, optional
            Id to denote new population. Default is None, and the id will be
            incremented from the previously simulated.
        
        Returns
        -------
        sc.AnnData
            Previously simulated single-cell dataset with new population
            appended.
        
        Raises
        ------
        ValueError
            Raised is provided `pop_id` is non-unique.
        """
        mus = adata.var['Base.Mu'].values
        if perturbed:
            mus *= adata.var['Perturbation.Shift'].values
        previous_markers = np.hstack([x for x in\
                                      population_markers(adata).values()])
        # don't set markers to perturbed genes, goal is to simulate an
        # unperturbed cell identity
        if perturbed:
            p_genes = adata.var.index[adata.var['Perturbation.Shift'] != 1]\
                           .values
            previous_markers = np.hstack([previous_markers, p_genes])
        n_markers = np.random.binomial(adata.shape[1],
                                       self.control_kwargs['p_marker'])
        new_markers = np.random.choice(list(set(range(adata.shape[1]))
                                            - set(previous_markers)),
                                       n_markers)
        pop_markers = np.array([False]*adata.shape[1])
        pop_markers[new_markers] = True
        gamma = stats.gamma(a=2, scale=2)
        mus[new_markers] *= gamma.rvs(len(new_markers))
        counts, __ = simulate_counts(n_cells, mus,
                                     adata.var['Base.Dispersion'].values,
                                     1, [n_cells])
        if pop_id is None:
            c_pops = adata.obs['Population'].unique()
            try:
                pop_id = c_pops.max() + 1
            except TypeError:
                pop_id = len(c_pops) + 1
        else:
            if pop_id in adata.obs['Population']:
                raise ValueError("Non-unique population id: {}".format(pop_id))
        pop_var = pd.DataFrame({'Pop.{}.Mu'.format(pop_id): mus,
                                'Pop.{}.Marker'.format(pop_id): pop_markers})
        pop_obs = pd.DataFrame({'Population': np.array([pop_id]*n_cells)})
        new_adata = sc.AnnData(X=counts, obs=pop_obs, var=pop_var)
        return utils.rbind_adata([adata, new_adata])
        
    def run(self, simulations=1, replications=1, controls=None,
            pop_targets=None):
        """
        Simulate control and perturbed datasets under experimental conditions.

        Parameters
        ----------
        simulations : int, optional
            Number of control datasets to simulate. Default is 1, and a single
            reference control dataset will be simulated.
        replications : int, optional
            Number of perturbations to simulate for each control dataset.
            Default is 1, and a single perturbation will be simulated for each
            reference control dataset.
        controls : sc.AnnData, optional
            A dataset of simulated control cells to perturb. Default is None,
            and a control dataset will be simulated according to parameters
            defined by `control_kwargs`.
        pop_targets : list-like, optional
            Populations to target during perturbation. Default is None, and
            perturbed genes will be randomly selected. If a list of populations
            is provided, marker genes for these populations will be targeted
            instead. 
        
        Returns
        -------
        list
            Two-dimensional list that is indexed first by simulation and second
            by replicate. 
        """
        out = []
        for __ in range(simulations):
            sim_out = []
            if controls is None:
                controls = self.simulate_controls()
            if not isinstance(controls, sc.AnnData):
                raise ValueError("Unexpected type for `controls`: {}".format(
                                  type(controls)))
            markers = population_markers(controls)
            if pop_targets is not None:
                self.perturb_kwargs['gene_targets'] = []
                for each in pop_targets:
                    self.perturb_kwargs['gene_targets'] += list(markers[each])
            else:
                markers = np.hstack([x for x in markers.values()])
                possible = set(controls.var.index).difference(markers)
                try:
                    p = self.perturb_kwargs['percent_perturb']
                except KeyError:
                    p = None
                if p is None:
                    p = 0.20
                targets = np.random.choice(list(possible),
                                           int(p * controls.shape[1]))
                self.perturb_kwargs['gene_targets'] = targets
                self.perturb_kwargs['percent_perturb'] = None
            for __ in range(replications):
                treated = perturb(controls, **self.perturb_kwargs)
                controls.obs['Treatment'] = 'Control'
                treated.obs['Treatment'] = 'Perturbed'
                combined = utils.rbind_adata([controls, treated])
                sim_out.append(combined)
            out.append(sim_out)
        return out


def average_exp(scale_factor, n=1):
    r"""
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
    r"""
    Estimate the probability of dropout for a given gene.

    Estimate the probability of a dropout even using a sigmoid function, and
    the following models. 

    ..math::
        p_i = 1 - sigmoid(x) \\
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
    return 1 - sigmoid(x)

def simulate_counts(n_samples, mus, dispersion, populations, pop_sizes,
                    percentile=50, method='conditional', dropout=0.66):
    r"""
    Simulate counts across genes for a set number of samples.
    
    Parameters
    ----------
    n_samples : int
        Number of samples to simulate.
    mus : np.ndarray
        A gene by population matrix, where each element represents the average
        expression value in the given population. 
    dispersion : numpy.ndarray
        A gene length vector with dispersion paramtersfor modelling counts
        using a negative binomial distribution for each gene.
    populations: int
        Number of populations to sample.
    pop_sizes : numpy.ndarray
        list of samples per population.
    percentile : float, optional
        Which percentile to use when calculating dropout probabilities. Default
        is 0.5, and the median will be used.
    method : str, optional
        How to perform dropout. Default is `conditional` and dropout rates for
        each gene will be conditioned on its average experession compared to the
        median average expression. The other option is `uniform` and all genes
        will share the same chance of dropping out. Default is conditional
    dropout : float, optional
        Value between 0 and 1 denoting the average chance of a dropout event for
        a given gene. Only used if `method=='uniform'`. Default is 0.66
    
    Returns
    -------
    (numpy.ndarray, numpy.ndarray)
        An :math:`n \times p` count matrix, where :math:`n` is number of
        samples, defined by `n_samples`, and :math:`p` is the number of genes,
        defined by the size of `mus`.

        A :math:`p` length vector of population labels for each row in returned
        count matrix.
    """
    # sample by gene expression matrix
    X_ = np.zeros((n_samples, mus.shape[0]))
    # vector labelling which population each sample belongs to
    labels_ = np.hstack([np.array([i + 1] * pop_sizes[i])\
                    for i in range(populations)])
    # ensure proper shape of dispersion vector
    if len(dispersion.shape) == 1:
        dispersion = dispersion.reshape((-1, 1))
    elif len(dispersion.shape) != 2\
    and dispersion.shape[0] < dispersion.shape[1]:
        raise ValueError("Expected column vector for dispersion. Received "
                         "vector with shape {}.".format(dispersion.shape))
    if len(mus.shape) == 1:
        mus = mus.reshape((-1, 1))
    # calculate theoretical expression averages across populations
    # |gene| x |populations| matrix
    means_ = mus \
           * np.ones((dispersion.shape[0], populations))\
           * dispersion
    percentile_ = np.percentile(means_, percentile)

    # calculate dropout probabilites for each gene in each population
    # a |gene| x |populations| size matrix
    if method == 'conditional':
        p_dropout = dropout_probability(means_, percentile_)
    elif method == 'uniform':
        p_dropout = np.ones_like(means_) * dropout
    # simulate counts across populations
    for i in range(populations):
        if i == 0:
            start = 0
        else:
            start = pop_sizes[:i].sum()
        # simulate gene counts for each cell in population i
        for j in range(mus.shape[0]):
            dist = stats.nbinom(dispersion[j, 0],
                                1 - mus[j, i] / (mus[j, i] + dispersion[j, 0]))
            # calculate probability of reads not dropping out with bernoulli
            # 1 = keep count, 0 = dropeed, multiplying by dropout vector
            # moves dropped counts to 0. 
            dropped = stats.bernoulli(p=1 - p_dropout[j, i]).rvs(pop_sizes[i])
            X_[start:start + pop_sizes[i], j] = dist.rvs(pop_sizes[i]) * dropped
    return X_, labels_
