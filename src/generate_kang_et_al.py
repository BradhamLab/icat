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

from downstream.src.analysis import utils as dutils

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
    idx_to_gene = {i + 1:x for i, x in enumerate(genes.index)}
    idx_to_bc = {i + 1:barcodes.loc[i, 'cell.barcode'] for i in barcodes.index}
    data = pd.read_csv(matrix_file, delimiter=' ', skiprows=3,
                       names=['gene', 'cell', 'count'])
    matrix = pd.pivot_table(data=data, index='cell', columns='gene',
                            values='count', fill_value=0)
    print(matrix.shape)
    matrix = matrix.rename(index=idx_to_bc, columns=idx_to_gene)
    return matrix


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        genes = pd.read_csv(snakemake.input['genes'],
                            names=['ensmbl.id', 'name'], delimiter='\t',
                            index_col=0)
        cells = pd.read_csv(snakemake.input['cells'],
                            index_col=0, delimiter='\t')
        mats_and_bars = [snakemake.input['mtx'], snakemake.input['barcodes']]
        count_mats = [create_count_matrix(mtx, bc, genes)\
                      for mtx, bc in zip(*mats_and_bars)]
        shared = set(count_mats[0].index.values)\
                .intersection(count_mats[1].index.values)
        count_mats[0] = count_mats[0].rename(index={x:x+'1' for x in shared})
        # combine count matrices, verify unique index ids, sort gene names. 
        counts = pd.concat(count_mats, axis=0, verify_integrity=True,
                           sort=True)
        counts.fillna(value=0, inplace=True)
        adata = sc.AnnData(X=counts.values, obs=cells.loc[counts.index, :],
                           var=genes.loc[counts.columns, :])
        adata = dutils.filter_cells(adata, 'multiplets',
                                    lambda x: x=='singlet').copy()
        adata.write_csvs(dirname=snakemake.params['outdir'], skip_data=False)

    
