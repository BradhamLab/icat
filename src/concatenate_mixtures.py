import pandas as pd
import numpy as np
from scanpy import api as sc
from icat.src import utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        X_csvs = sorted(snakemake.input['X'])
        obs_csvs = sorted(snakemake.input['obs'])
        var_csvs = sorted(snakemake.input['var'])
        adatas = []
        for X, obs, var in zip(X_csvs, obs_csvs, var_csvs):
            X = np.loadtxt(X, delimiter=',')
            obs = pd.read_csv(obs, index_col=0)
            var = pd.read_csv(var, index_col=0)
            adatas.append(sc.AnnData(X=X, obs=obs, var=var))
        combined = utils.rbind_adata(adatas)
        combined.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)