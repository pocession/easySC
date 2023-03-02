## library

import argparse
from pathlib import Path

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

    p.data = p.data.resolve() # full file path
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

    return paths

def load_data(paths):
    '''
    Load the 3 data files
    Save to a data class object
    '''
    print(paths)

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
