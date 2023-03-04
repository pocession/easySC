## library

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc


## Loading packages
print("..Loading packages..")
'''
scanpy read all three files and put it in an self.data object. Not sure if we still need a class object here.
self.data: https://anndata.readthedocs.io/en/latest/
scanpy: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.htm
'''
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

class SCAnalysis(paths):
    def __init__(paths):
        self.paths = paths

    def load_data(self):
        '''Load the 3 data files'''

        print("..Reading input data..")
        self.data = sc.read_10x_mtx(
        self.paths[3],  # the directory with the `.mtx` file
        var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
        cache=True)
        self.data.var_names_make_unique() # necessary because we use gene symbols as var_names
        print(f"Your data set contains {self.data.n_obs} cells and {self.data.n_vars} genes.")

    def filter_data(self):
        '''ToDo: perform basic filtering, umi/gene, gene/cell'''
        
        # basic filtering
        ## We should make the app can take the criteria from input
        sc.pp.filter_cells(self.data, min_genes=200)
        sc.pp.filter_genes(self.data, min_cells=3)
        print(f"After filtering, your data set contains {self.data.n_obs} cells and {self.data.n_vars} genes.")

    def qc_data(self):
        ''' 
        Draw figures for raw data and save them.
        They are useful for QC.
        The draw function comes from scanpy itself.
        The out plots will be save in "./figures" automatically.
        '''

        self.data.var['mt'] = self.data.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(self.data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(self.data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                        jitter=0.4, multi_panel=True, show=False, save="_basicQC.pdf")
        sc.pl.scatter(self.data, x='total_counts', y='pct_counts_mt',show=False, save="_mtQC.png")
        sc.pl.scatter(self.data, x='total_counts', y='n_genes_by_counts',show=False, save="_geneQC.png")

    def make_analysis(self):
        '''generate analysis and plot'''
        pass

    def save_csv(self):
        '''ToDo: write results into csv'''
        pass

## main
if __name__ == "__main__":

    paths = get_input()
    exp = SCAnalysis(paths)
    exp.load_data()
    exp.filter_data()
    exp.qc_data()
    exp.make_analysis()
    exp.save_csv()