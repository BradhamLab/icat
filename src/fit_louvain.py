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
    data_dir = '../data/processed'
