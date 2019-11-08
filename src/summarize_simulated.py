import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from icat.src import plotting
from cycler import cycler

import os

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        if not os.path.exists(snakemake.params['plotdir']):
            os.makedirs(snakemake.params['plotdir'])
        results = pd.read_csv(snakemake.input['results'], index_col=0)
        metrics = list(plotting.metric_dictionary.keys())
        results[['Experiment', 'Perturbation', 'Sim', 'Rep']] =\
            results[['Experiment', 'Perturbation', 'Sim', 'Rep']].fillna('')
        # create column to group runs 
        results['ExpID'] = results[['Experiment', 'Perturbation',
                                    'Sim', 'Rep']].apply(lambda x: "".join(x),
                                                         axis=1)
        # set index to always keep track of methods and runs
        results.index = pd.MultiIndex.from_arrays([results.index,
                                                   results['method'],
                                                   results['ExpID']])
        results.index.names = [None, 'Method', 'Run']
        # rank performance in metrics between methods across simulations and
        # replications
        ranks = results.groupby(['Run'])[metrics].rank(ascending=False)
        ranks = ranks.groupby('Method').mean().rank()
        # plot average of average ranks
        plotting.ranked_heatmap(ranks)
        plt.savefig(snakemake.output['ranks'])
        plotting.close_plot()
        # average performance across replications
        results['ExpID'] = results[['Experiment', 'Perturbation', 'Sim']].apply(
                                    lambda x: "".join(x), axis=1)
        sim_mean = results.groupby(['ExpID', 'method']).mean()
        # get Experiment averages across simulations
        sim_mean['id'] = sim_mean.apply(lambda x: x.name[0].split('Sim')[0],
                                        axis=1)
        by_exp = sim_mean.groupby(['id', 'method'])
        method_mean = by_exp.mean()
        method_deviation = by_exp.std()
        # fill na with 0's b/c not sure what happens lmao
        method_deviation.fillna(0, inplace=True)
        # plot performances across experiments
        plt.rcParams['axes.prop_cycle'] = cycler(color=plotting.method_colors)
        for exp in method_mean.index.get_level_values('id').unique():
            write_dir = os.path.join(snakemake.params['outdir'], exp)
            method_mean.loc[exp].to_csv(os.path.join(write_dir,
                                                     'metric_means.csv'))
            method_deviation.loc[exp].to_csv(os.path.join(write_dir,
                                                          'metric_stds.csv'))
            plotting.metric_plot(method_mean.loc[exp].T,
                                 method_deviation.loc[exp].T) 
            plt.savefig(os.path.join(snakemake.params['plotdir'],
                        exp + '_metrics.svg'))
            plotting.close_plot()
        