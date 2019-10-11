import pandas as pd
import numpy as np

import os

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
        sim_mean = results.groupby(['ExpID', 'method']).mean()
        sim_mean['id'] = sim_mean.apply(lambda x: x.name[0].split('Sim')[0],
                                        axis=1)
        by_exp = sim_mean.groupby(['id', 'method'])
        method_mean = by_exp.mean()
        method_deviation = by_exp.std()
        # fill na with 0's b/c not sure what happens lmao
        method_deviation.fillna(0, inplace=True)
        for exp in method_mean['id']:
            method_mean.loc[exp].to_csv(os.path.join(
                                            snakemake.params['outdir']),
                                            exp + '_scores.csv')
            method_deviation.loc[exp].to_csv(os.path.join(
                                                 snakemake.params['outdir']),
                                                 exp + '_devs.csv')
        