import pandas as pd

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        files = snakemake.input['csvs']
        dfs = []
        for each in files:
            dfs.append(pd.read_csv(each, index_col=0))
        out_df = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True)
        out_df.to_csv(snakemake.output['csv'])
