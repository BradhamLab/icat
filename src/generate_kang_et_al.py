"""
Functions to create single-cell datasets.

Combines data from Kang et al. 2018.

author: Dakota Y. Hawkins
contact: dyh0110@bu.edu
"""
import pandas as pd
from scanpy import api as sc
from glob import glob
import os
import pickle as pkl

def match_matrix_and_barcode_files(data_dir):
    """
    Match cell expression data files with associated cell barcode files.
    
    Parameters
    ----------
    data_dir : str
        Directory containing cell expression data and barcode files. Expression
        data should have *.mtx extensions, while barcode files are expected to
        have *.tsv extensions.
    
    Returns
    -------
    list, tuple
        List of tuples where each entry is a matched pair of expression and
        barcode files. Each element is ordered with the expression file first, 
        and barcode file second.
    
    Raises
    ------
    IOError
        Raised if a matrix file cannot be matched with a barcode file.
    """
    matrix_files = glob(os.path.join(data_dir, '*.mtx'))
    ids = [os.path.basename(x.split('_')[0]) for x in matrix_files]
    barcodes = [data_dir + x + '_barcodes.tsv' for x in ids]
    for i, each in enumerate(barcodes):
        if not os.path.exists(each):
            raise IOError("Cannot find matching barcode file "
                          "for {}. Expected {}".format(matrix_files[i], each))
    return [(matrix_files[i], barcodes[i]) for i in range(len(matrix_files))]


def create_count_matrix(matrix_file, barcode_file, genes):
    """
    Create a count matrix from necessary files.
    
    Parameters
    ----------
    matrix_file : str
        A sparse matrix file with 3 columns separated by spaces. In order, the
        expected columns are gene number, cell number, and umi count.
    barcode_file : str
        A single-column file containing all cell barcodes. Expected to be in
        order, such that "cell 1" listed in the matrix file will be the first
        entry.
    genes : pd.DataFrame
        A pandas dataframe containing gene annotations. Order is expected, such
        that 'gene 1' in the matrix file will be the first indexed gene.
    
    Returns
    -------
    pd.DataFrame
        An (n x p) data frame where n is the number of cells and p is the
        number of genes with at least one non-zero read.
    """
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


def main(data_dir, outdir):
    """
    Create an annotated dataframe object from the Kang 2018 dataset.
    
    Parameters
    ----------
    data_dir : str
        Directory containing downloaded data.
    outfile : str
        File to write pickled AnnData object to.
    """
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
    anno_df.write_csvs(dirname=outdir, skip_data=False)


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        main(snakemake.params['datadir'], snakemake.params['outdir'])


    
