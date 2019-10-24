import numpy as np
import pandas as pd
import os

from icat.src import utils

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        csv_dict = {}
        for csv in snakemake.input['Xs']:
            run = os.path.basename(os.path.split(csv)[0])
            X = np.loadtxt(csv, delimiter=',')
            dropout = np.sum(X == 0) / X.size 
            exp_data = utils.parse_sim(run)
            exp_data['dropout'] = dropout
            csv_dict[run] = exp_data
        df = pd.DataFrame(csv_dict).T
        df.to_csv(snakemake.output['csv'])
