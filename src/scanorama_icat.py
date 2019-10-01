import pandas as pd
import numpy as np
from scanpy import api as sc
from downstream.src.analysis import utils as dutils
from icat.src import utils
from icat.src import models
import json

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        integrated = sc.AnnData(X=np.loadtxt(snakemake.input['X'],
                                             delimiter=','),
                                obs=pd.read_csv(snakemake.input['obs'],
                                                index_col=0),
                                var=pd.read_csv(snakemake.input['var'],
                                                index_col=0))
        integrated.var.index = integrated.var.index.astype(str)
        integrated.obs.index = integrated.obs.index.astype(str)
        isolated = dutils.filter_cells(integrated, snakemake.params['treat_col'],
                                      lambda x: x not in\
                                            snakemake.params['treat_values']).\
                                            copy()
        mixed = dutils.filter_cells(integrated, snakemake.params['treat_col'],
                                      lambda x: x in\
                                             snakemake.params['treat_values']).\
                                             copy()
        with open(snakemake.input['ncfs'], 'r') as f:
            icat_kws = json.load(f) 
        with open(snakemake.input['json'], 'r') as f:
            fit_data = json.load(f)
        icat_kws['neighbor_kws']['n_neighbors'] = fit_data['n_neighbors']
        icat_model = models.icat(**icat_kws)
        out = icat_model.cluster(isolated, mixed)
        out.obs.rename(columns={'sslouvain': 'scanorama.sslouvain'},
                       inplace=True)
        out.write_csvs(dirname=snakemake.params['outdir'],
                       skip_data=False)
        
        