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
        prtb_adatas = []
        # separate control dataset from mixed cells
        for i, each in enumerate(all_X):
            if snakemake.params['control_id'] in each:
                # load in control data
                ctrl_adata = sc.AnnData(X=np.loadtxt(each, delimiter=','),
                                        obs=pd.read_csv(all_obs[i], index_col=0))
            else:
                # load in mixed data
                prtb_adatas.append(sc.AnnData(X=np.loadtxt(each, delimiter=','),
                                              obs=pd.read_csv(all_obs[i],
                                                              index_col=0)))
        # combine adatas
        prtb_combined = utils.rbind_adata(prtb_adatas)
        with open(snakemake.input['ncfs'], 'r') as f:
            icat_kws = json.load(f)
        with open(fit_json, 'r') as f:
            fit_data = json.load(f)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
        # cluster data using icat
        icat_model = models.icat(**icat_kws)
        combined = icat_model.cluster(ctrl_adata, prtb_combined)
        combined.write_csvs(dirname=snakemake.params['outdir'],
                            skip_data=False)
        