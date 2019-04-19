import pandas as pd
from scanpy import api as sc
from glob import glob
import os
import pickle as pkl

def match_matrix_and_barcode_files(data_dir):

    matrix_files = glob(os.path.join(data_dir, '*.mtx'))
    ids = [os.path.basename(x.split('_')[0]) for x in matrix_files]
    barcodes = [data_dir + x + '_barcodes.tsv' for x in ids]
    for i, each in enumerate(barcodes):
        if not os.path.exists(each):
            raise IOError("Cannot find matching barcode file "
                          "for {}. Expected {}".format(matrix_files[i], each))
    return [(matrix_files[i], barcodes[i]) for i in range(len(matrix_files))]


def create_count_matrix(matrix_file, barcode_file, genes):
    # read in cell barcodes
    barcodes = pd.read_csv(barcode_file, names=['cell.barcode'])
    # match index numbers to names, indices started with 1
    idx_to_gene = {i + 1:x for i, x in
                   enumerate(genes.index)}
    idx_to_bc = {i + 1:barcodes.loc[i, 'cell.barcode'] for i in barcodes.index}
    data = pd.read_csv(matrix_file, delimiter=' ', skiprows=3,
                       names=['gene', 'cell', 'count'])
    matrix = pd.pivot_table(data=data, index='cell', columns='gene',
                            values='count', fill_value=0)
    matrix = matrix.rename(index=idx_to_bc, columns=idx_to_gene)
    return matrix


def main(data_dir, outfile):
    if data_dir[-1] != '/':
        data_dir += '/'
    files = match_matrix_and_barcode_files(data_dir)
    gene_file = glob(data_dir + '*genes.tsv')[0]
    cell_file = glob(data_dir + '*tsne*tsv')[0]
    genes = pd.read_csv(gene_file, names=['ensmbl.id', 'name'], delimiter='\t',
                        index_col=0)
    cells = pd.read_csv(cell_file, index_col=0, delimiter='\t')
    count_mats = [create_count_matrix(each[0], each[1], genes)\
                  for each in files]
    shared = set(count_mats[0].index.values)\
             .intersection(count_mats[1].index.values)
    count_mats[0] = count_mats[0].rename(index={x:x+'1' for x in shared})
    # combine count matrices, verify unique index ids, sort gene names. 
    counts = pd.concat(count_mats, axis=0, verify_integrity=True,
                       sort=True)
    counts.fillna(value=0)
    anno_df = sc.AnnData(X=counts, obs=cells,
                         var=genes.loc[counts.columns.values, :])
    with open(outfile, 'wb') as f:
        pkl.dump(anno_df, f)


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        main(snakemake.params[0], snakemake.output[0])


    
