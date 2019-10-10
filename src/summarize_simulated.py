import pandas as pd
import numpy as np

import seaborn as sns

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        results = pd.read_csv(snakemake.input['results'], index_col=0)
        results[['Experiment', 'Perturbation', 'Sim', 'Rep']] =\
            results[['Experiment', 'Perturbation', 'Sim', 'Rep']].fillna('')
        results['ExpID'] = results.apply(lambda x: "{}{}{}".format(
                                            x['Experiment'],
                                            x['Perturbation'],
                                            x['Sim']),
                                         axis=1)
        method_mean = results.groupby(['ExpID', 'method']).mean()
        method_deviation = results.groupby(['ExpID', 'method']).std()

        metrics = ['adjusted.mutual.info', 'homogeneity', 'completeness',
                   'adjusted.rand', 'fowlkes.mallows']
        g = sns.catplot(x='method', y=metrics, )