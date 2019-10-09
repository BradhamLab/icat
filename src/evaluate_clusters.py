import pandas as pd
from icat.src import utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        method_to_label = {x: y for x, y in zip(['icat', 'seurat',
                                                 'scanorama', 'icat_scan'],
                                                ['NCFS-SSLouvain', 'cluster',
                                                 'scanorama.louvain',
                                                 'scanorama.sslouvain'])}
        labels = snakemake.params['identity']
        method = snakemake.params['method']
        data = pd.read_csv(snakemake.input['obs'], index_col=0)
        performance = utils.performance(data, labels, method_to_label[method])
        performance['method'] = method
        exp_info = utils.parse_sim(snakemake.params['exp'])
        for key, value in exp_info.items():
            performance[key] = value
        out = pd.DataFrame(performance, index=[snakemake.params['exp']])
        out.to_csv(snakemake.output['csv'])