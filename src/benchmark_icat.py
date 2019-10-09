import pandas as pd
import numpy as np
from scanpy import api as sc
from icat.src import models, utils
from downstream.src.analysis import utils as dutils
import json

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        fit_json = snakemake.input['json']
        adata = sc.AnnData(X=np.loadtxt(snakemake.input['X'], delimiter=','),
                           obs=pd.read_csv(snakemake.input['obs'], index_col=0),
                           var=pd.read_csv(snakemake.input['var'], index_col=0))
        ctrls = dutils.filter_cells(adata, snakemake.params['treatment'],
                                    lambda x: x == snakemake.params['control'])
        prtbs = dutils.filter_cells(adata, snakemake.params['treatment'],
                                    lambda x: x != snakemake.params['control'])
        # combine adatas
        with open(snakemake.input['ncfs'], 'r') as f:
            icat_kws = json.load(f)
        with open(fit_json, 'r') as f:
            fit_data = json.load(f)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
        # cluster data using icat
        icat_model = models.icat(**icat_kws)
        combined = icat_model.cluster(ctrls, prtbs)
        combined.write_csvs(dirname=snakemake.params['outdir'],
                            skip_data=False)
        