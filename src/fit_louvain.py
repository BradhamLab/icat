from sklearn import metrics
import numpy as np
import scanpy.api as sc

def main(andata, label_col):
    sc.pp.pca(andata)
    performance = {'score': -np.inf,
                   'n_neighbors': None,
                   'resolution': None}
    n_max = int(0.75 * andata.shape[0])
    n_min = 3
    resolutions = [1]
    for n in range(n_min, n_max):
        for r in resolutions:
            sc.pp.neighbors(andata, n_neighbors=n)
            sc.tl.umap(andata, min_dist=0.0)
            sc.tl.louvain(andata, resolution=r, key_added='louvain')
            score = metrics.adjusted_rand_score(andata.obs[label_col],
                                                andata.obs['louvain'])
            if score > performance['score']:
                performance['score'] = score
                performance['n_neighbors'] = n
                performance['resolution'] = r
    return performance

if __name__ == '__main__':
    from glob import glob
    import pandas as pd
    import pickle as pkl
    import os
    data_dir = '../data/processed'
    label_col = 'Population'
    out_csv = '../data/interim/control_louvain_fits.csv'
    regex = '*Controls.pkl'
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        data_dir = snakemake.params['data']
        label_col = snakemake.params['label']
        out_csv = snakemake.output['csv']
        regex = snakemake.params['regex']
    datasets = sorted(glob(os.path.join(data_dir, regex)))
    performance = dict()
    for fn in datasets:
        with open(fn, 'rb') as f:
            adata = pkl.load(f)
        name = os.path.basename(fn)
        performance[name] = main(adata, label_col)
    pd.DataFrame(performance).T.to_csv(out_csv)