import pandas as pd
import numpy as np
from scanpy import api as sc
from icat.src import models, utils
import json



if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        fit_json = snakemake.input['json']
        all_X = sorted(snakemake.input['X'])
        all_obs = sorted(snakemake.input['obs'])
        ctrl_X = None
        ctrl_obs = None
        for i, each in enumerate(all_X):
            if snakemake.params['control_id'] in each:
                ctrl_X = each
                ctrl_obs = all_obs[i]
                all_X.remove(ctrl_X)
                all_obs.remove(ctrl_obs)
        ctrl_adata = sc.AnnData(X=np.loadtxt(ctrl_X, delimiter=','),
                                obs=pd.read_csv(ctrl_obs, index_col=0))
        prtb_adatas = []
        for i in range(len(all_X)):
            prtb_adatas.append(sc.AnnData(X=np.loadtxt(all_X[i],
                                                       delimiter=','),
                                          obs=pd.read_csv(all_obs[i],
                                                          index_col=0)))
        prtb_combined = utils.rbind_adata(prtb_adatas)
        # 
        icat_kws = {'ncfs_kws': {'sigma': 2, 'reg': 1,
                                 'kernel': 'exponential'},
                    'neighbor_kws': {'n_neighbors': None},
                    'cluster_kws': {'resolution': None},
                    'weight_threshold': 1.0,
                    'treatment_col': 'benchmark'}
        with open(fit_json, 'r') as f:
            fit_data = json.load(f)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
        icat_model = models.icat(**icat_kws)
        combined = icat_model.cluster(ctrl_adata, prtb_combined)
        combined.write_csvs(dirname=snakemake.params['outdir'],
                            skip_data=False)
        