## library

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

print("..Loading packages..")
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
print("..Loading packages complete!")

## functions
def get_input():
    '''
    Read data and check input
    https://docs.python.org/3/library/argparse.html
    '''
    parser = argparse.ArgumentParser(description='Generate a basic profiles of your single cell RNA-seq data.')
    parser.add_argument('--data','-d', 
                        required=True,
                        type=Path,
                        help='specify the input path') 
    p = parser.parse_args()

    # p.data = input folder
    # p.data.resolve() = full path of input folder
    p.data = p.data.resolve()
    if p.data.is_dir():
        print(f"Input data dir: {p.data}")
    else:
        print(f"Does not exist: {p.data}")
    
    files = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]
    paths = list()
    for f in files:
        f = p.data / f
        if f.is_file():
            paths.append(f)
        else:
            print(f"missing file: {f}")
    paths.append(p.data) # append the input folder path in the last position

    return paths

def load_data(paths):
    '''
    Load the 3 data files
    Save to a data class object

    scanpy read all three files and put it in an adata object. Not sure if we still need a class object here.
    adata: https://anndata.readthedocs.io/en/latest/
    scanpy: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
    '''

    # print(paths)
    print("..Reading input data..")
    adata = sc.read_10x_mtx(
    paths[3],  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
    adata.var_names_make_unique() # necessary because we use gene symbols as var_names
    print(f"Your data set contains {adata.n_obs} cells and {adata.n_vars} genes.")
    return adata

def filter_data(adata):
    '''ToDo: perform basic filtering, umi/gene, gene/cell'''
    
    # basic filtering
    ## We should make the app can take the criteria from input
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After filtering, your data set contains {adata.n_obs} cells and {adata.n_vars} genes.")
    return adata

def qc_data(adata_filtered):
    ''' 
    Draw figures for raw data and save them.
    They are useful for QC.
    The draw function comes from scanpy itself.
    The out plots will be save in "./figures" automatically.

    '''
    adata = adata_filtered

    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                    jitter=0.4, multi_panel=True, show=False, save="_basicQC.pdf")
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',show=False, save="_mtQC.png")
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',show=False, save="_geneQC.png")

def make_analysis(adata):
    '''generate analysis and plot'''
    pass

def save_csv():
    '''ToDo: write results into csv'''
    pass

## main
if __name__ == "__main__":

    paths = get_input()
    adata = load_data(paths)
    adata_filtered = filter_data(adata)
    qc_data(adata_filtered)
    make_analysis(adata)

    save_csv()