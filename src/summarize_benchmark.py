import pandas as pd
from icat.src import utils
from matplotlib import pyplot as plt
from icat.src import plot_performance
import re
import os

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        performances = {}
        prefix = os.path.commonpath(snakemake.input['obs']).replace('/', '\/')
        method_regex = re.compile('(?<={}\/)(.*?)(?=\/)'.format(prefix))
        for each in snakemake.input['obs']:
            method = method_regex.search(each).group()
            data = pd.read_csv(each, index_col=0)
            print(method)
            performances[method] = utils.performance(data,
                                               snakemake.params['identity'],
                                               utils.label_dictionary[method])
        out = pd.DataFrame(performances)
        out.to_csv(snakemake.output['csv'])

        plot_performance.metric_plot(out)
        plt.savefig(snakemake.output['svg'])
