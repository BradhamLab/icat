import pandas as pd
import scanpy.api as sc
from icat.src import utils
import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        clustered = pd.read_csv(snakemake.input['csv'], index_col=0)
        X = np.zeros((clustered.shape[0], 1))
        adata = sc.AnnData(X=X, obs=clustered)
        adata.obsm['X_umap'] = adata.obs[['UMAP1', 'UMAP2']].values
        perf = utils.performance(adata, 'Population', 'cluster')
        perf['method'] = 'seurat'
        sim_data = utils.parse_sim(snakemake.params['name'])
        out_df = pd.DataFrame(dict(perf, **sim_data),
                              index=[snakemake.params['name']])
        out_df.to_csv(snakemake.output['csv'])
        utils.plot_umap(adata, 'cluster', 'Population')
        plt.savefig(snakemake.output['svg'])
