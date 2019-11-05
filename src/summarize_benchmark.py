import os
import re

import pandas as pd
from matplotlib import pyplot as plt
from cycler import cycler

from icat.src import plot_performance, utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        identity = snakemake.params['identity']
        plotdir = snakemake.params['plotdir']
        treatment = snakemake.params['treatment']
        control_id = snakemake.params['controls']
        performances = {}
        prefix = os.path.commonpath(snakemake.input['obs']).replace('/', r'\/')
        method_regex = re.compile(r'(?<={}\/)(.*?)(?=\/)'.format(prefix))
        plotted_expectation = False
        for each in snakemake.input['obs']:
            method = method_regex.search(each).group()
            data = pd.read_csv(each, index_col=0)
            column = utils.label_dictionary[method]
            performances[method] = utils.performance(data, identity, column)
            # subset data to control and perturbed to asses performance within
            # and out of treatments
            ctrls = data[data[treatment] == control_id]
            prtbs = data[data[treatment] != control_id]
            for x, subset in zip(['all', 'controls', 'treated'],
                                 [data, ctrls, prtbs]):
                # ternary plots
                plot_performance.ternary_plot(data, column, identity)
                plt.savefig(os.path.join(plotdir, '{}_{}_ternary.svg'.\
                                                  format(x, method)))
                plot_performance.close_plot()
                if not plotted_expectation:
                    plot_performance.ternary_plot(data, identity, identity)
                    plt.savefig(os.path.join(plotdir,
                                             '{}_expected_ternary.svg'.\
                                             format(x)))
            plotted_expectation = True
        out = pd.DataFrame(performances)
        out.to_csv(snakemake.output['csv'])
        plt.rcParams['axes.prop_cycle'] =\
            cycler(color=plot_performance.method_colors)
        plot_performance.metric_plot(out)
        plt.savefig(snakemake.output['svg'])
