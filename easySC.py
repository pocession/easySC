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
    '''
    print(paths)
    adata = sc.read_10x_mtx(
    paths[3],  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
    adata

def filter_data():
    '''ToDo: perform basic filtering, umi/gene, gene/cell'''
    pass

def make_analysis():
    '''generate analysis and plot'''
    pass

def save_csv():
    '''ToDo: write results into csv'''
    pass

def save_plot():
    '''save plots'''
    pass


## main
if __name__ == "__main__":

    paths = get_input()

    load_data(paths)
    filter_data()
    make_analysis()

    save_csv()
    save_plot()
