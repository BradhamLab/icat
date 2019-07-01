import pickle as pkl
import scanpy.api as sc

def read_write(pkl_file, outdir):
    with open(pkl_file, 'rb') as f:
        data = pkl.load(f)
    data.write_csvs(outdir, skip_data=False)

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        read_write(snakemake.input['ctrl'], snakemake.params['ctrl_dir'])
        read_write(snakemake.input['prtb'], snakemake.params['prtb_dir'])