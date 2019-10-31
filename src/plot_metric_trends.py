import os
import re

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns

from icat.src import plot_performance

def experiment_match(exp_id, index):
    regex = re.compile(exp_id)
    return any([regex.search(x) is not None for x in index])

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        results = pd.read_csv(snakemake.input['results'], index_col=0)
        simulations = pd.read_csv(snakemake.input['simulations'], index_col=0)
        exp_id = snakemake.params['exp']
        feature = snakemake.params['sim_feature']
        subset = simulations[simulations.apply(lambda x: exp_id in x.name,
                                               axis=1)]
        if feature == 'dropout':
            results['ExpID'] = results.apply(lambda x: '{}.{}{}{}'.format(
                                                                x.Experiment,
                                                                x.Perturbation,
                                                                x.Sim,
                                                                x.Rep),
                                             axis=1)
        else:
            # match ExpID column in results Dataframe to subset index values
            results['ExpID'] = results.apply(lambda x: '{}.{}'.format(x.Experiment,
                                                                    x.Perturbation),
                                            axis=1)
        results = results[results.apply(lambda x: x.ExpID in subset.index,
                                        axis=1)]
        results.loc[results.index, feature] = results.apply(lambda x:\
                                                   subset.loc[x.ExpID, feature],
                                                   axis=1)
        results.loc[results.index, 'ExpID'] = results.apply(lambda x: "{}{}{}"\
                                                     .format(x['Experiment'],
                                                             x['Perturbation'],
                                                             x['Sim']),
                                                            axis=1)
        sim_mean = results.groupby(['ExpID', 'method']).mean()
        sim_mean['Method'] = sim_mean.index.get_level_values('method')
        sim_mean['Method'] = sim_mean['Method'].apply(lambda x:\
                                          plot_performance.method_dictionary[x])
        sim_mean.rename(columns=plot_performance.metric_dictionary,
                        inplace=True)
        to_metric = {v:k for k, v in plot_performance.metric_dictionary.items()}
        sim_mean[feature] = sim_mean[feature] * 100
        for y in plot_performance.metric_dictionary.values():
            plot_performance.trendplot(sim_mean, x=feature,
                                       y=y, hue='Method',
                                       xlabel=snakemake.params['xlabel'])
            plt.savefig(os.path.join(snakemake.params['plotdir'],
                        to_metric[y].replace('.', '-') + '.svg'))
            plot_performance.close_plot()
        sim_mean.to_csv(snakemake.output['csv'])
