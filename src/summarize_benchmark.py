import pandas as pd
from icat.src import utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        keys = ['icat', 'seurat', 'scanorama', 'icat_scan']
        columns = ['sslouvain', 'cluster', 'scanorama.louvain',
                   'scanorama.sslouvain']
        performances = {}
        for i, each in enumerate(keys):
            data = pd.read_csv(snakemake.input[each], index_col=0)
            # print(data.columns)
            performances[each] = utils.performance(data,
                                                   snakemake.params['identity'],
                                                   columns[i])
        out = pd.DataFrame(performances).T
        out.to_csv(snakemake.output['csv'])