from sklearn import metrics
import numpy as np
import pandas as pd
import scanpy.api as sc
import matplotlib.pyplot as plt
import pickle as pkl
import os

def main(andata, label_col):
    sc.pp.pca(andata)
    performance = {'score': -np.inf,
                   'n_neighbors': None,
                   'resolution': None}
    n_max = int(0.60 * andata.shape[0])
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

def plot_cells(andata, n, fn):
    sc.pp.pca(andata)
    sc.pp.neighbors(andata, n_neighbors=n)
    sc.tl.umap(andata, min_dist=0)
    sc.pl.umap(andata, color='Population', save=fn, show=False)

if __name__ == '__main__':
    from glob import glob
    import pandas as pd
    import pickle as pkl
    import os
    data_dir = '../data/processed'
    label_col = 'Population'
    out_csv = '../data/interim/control_louvain_fits.csv'
    regex = '*Controls.pkl'
    plot_dir = '../plots/simulated'
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        data_dir = snakemake.params['data']
        label_col = snakemake.params['label']
        out_csv = snakemake.output['csv']
        regex = snakemake.params['regex']
        plot_dir = snakemake.params['plotdir']
    sc.settings.figdir = plot_dir
    datasets = sorted(glob(os.path.join(data_dir, regex)))
    performance = dict()
    for fn in datasets:
        with open(fn, 'rb') as f:
            adata = pkl.load(f)
        with open(fn.replace('Controls', 'Treated'), 'rb') as f:
            perturbed = pkl.load(f)
        name = os.path.basename(fn)
        performance[name] = main(adata, label_col)
        combined = sc.AnnData(X=np.vstack((adata.X, perturbed.X)),
                              obs=pd.concat((adata.obs, perturbed.obs)),
                              var=adata.var)
        base = name.replace('Controls.pkl', "")
        names = [base + x for x in ['Controls.png', 'Perturbed.png',
                                    "Combined.png"]]
        for data, plotfile in zip([adata, perturbed, combined], names):
            plot_cells(data, int(performance[name]['n_neighbors']), plotfile)
            plt.cla()
            plt.clf()
            plt.close()
    pd.DataFrame(performance).T.to_csv(out_csv)