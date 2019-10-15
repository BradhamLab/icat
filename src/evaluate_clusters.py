import pandas as pd
from icat.src import utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        method_to_label = {x: y for x, y in zip(['icat', 'seurat233',
                                                 'seurat311',
                                                 'scanorama', 'icat_scan',
                                                 'seurat_icat'],
                                                ['NCFS-SSLouvain', 'cluster',
                                                 'seurat_clusters',
                                                 'scanorama.louvain',
                                                 'scanorama.sslouvain',
                                                 'seurat.sslouvain'])}
        labels = snakemake.params['identity']
        method = snakemake.params['method']
        data = pd.read_csv(snakemake.input['obs'], index_col=0)
        performance = utils.performance(data, labels, method_to_label[method])
        performance['method'] = method
        run_info = utils.parse_sim(snakemake.params['run'])
        for key, value in run_info.items():
            performance[key] = value
        out = pd.DataFrame(performance, index=[snakemake.params['run']])
        out.to_csv(snakemake.output['csv'])