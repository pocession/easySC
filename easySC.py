## library

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

print("..Loading packages..")
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
print("..Loading packages complete!")

## variables
# This variables should be taken from input also
min_genes = 200
min_cells = 3
up_cutoff = 0.75
down_cutoff = 0.25

## functions
def get_input():
    """
    Read data and check input
    https://docs.python.org/3/library/argparse.html
    """
    parser = argparse.ArgumentParser(
        description="Generate a basic profiles of your single cell RNA-seq data."
    )
    parser.add_argument(
        "--data", "-d", required=True, type=Path, help="specify the input path"
    )
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
    paths.append(p.data)  # append the input folder path in the last position

    return paths


class SCAnalysis:
    def __init__(self, paths):
        self.paths = paths

    @classmethod
    def load_data(exp, paths):
        """
        Load the 3 data files
        Save to a data class object

        scanpy read all three files and put it in an adata object. Not sure if we still need a class object here.
        adata: https://anndata.readthedocs.io/en/latest/
        scanpy: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
        """

        # print(paths)
        print("..Reading input data..")
        adata = sc.read_10x_mtx(
            paths[3],  # the directory with the `.mtx` file
            var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
            cache=True,
        )
        adata.var_names_make_unique()  # necessary because we use gene symbols as var_names
        print(f"Your data set contains {adata.n_obs} cells and {adata.n_vars} genes.")
        return adata

    def filter_data(exp, adata):
        """

        ToDo: perform filtering and produce QC plots

        """

        # Preliminary filtering
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(
            f"Perform the preliminary filtering process. Cells with less than {min_genes} expressed genes and genes not expressed in more than {min_cells} are excluded."
        )
        print(
            f"After the preliminary filtering process, your data set contains {adata.n_obs} cells and {adata.n_vars} genes."
        )

        # Draw figures for raw data and save them.
        # The draw function comes from scanpy itself.
        # # The out plots will be save in "./figures" automatically.

        adata.var["mt"] = adata.var_names.str.startswith(
            "MT-"
        )  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        ## This should be wrapped up into a plot function
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False,
            save="_preliminaryQC.pdf",
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="pct_counts_mt",
            show=False,
            save="_preliminaryQC_mt.png",
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="n_genes_by_counts",
            show=False,
            save="_preliminaryQC_gene.png",
        )

        # Determine the cutoff limits
        df_obs = pd.DataFrame(adata.obs)
        pct_counts_mt_up_cutoff = df_obs["pct_counts_mt"].quantile(up_cutoff)
        total_counts_up_cutoff = df_obs["total_counts"].quantile(up_cutoff)
        n_genes_by_count_up_cutoff = df_obs["n_genes_by_counts"].quantile(up_cutoff)

        pct_counts_mt_down_cutoff = df_obs["pct_counts_mt"].quantile(down_cutoff)
        total_counts_down_cutoff = df_obs["total_counts"].quantile(down_cutoff)
        n_genes_by_count_down_cutoff = df_obs["n_genes_by_counts"].quantile(down_cutoff)

        # Further filtering process
        adata = adata[adata.obs.total_counts < total_counts_up_cutoff, :]
        adata = adata[adata.obs.n_genes_by_counts < n_genes_by_count_up_cutoff, :]
        adata = adata[adata.obs.pct_counts_mt < pct_counts_mt_up_cutoff, :]

        adata = adata[adata.obs.total_counts > total_counts_down_cutoff, :]
        adata = adata[adata.obs.n_genes_by_counts > n_genes_by_count_down_cutoff, :]
        adata = adata[adata.obs.pct_counts_mt > pct_counts_mt_down_cutoff, :]

        ## This should be wrapped up into a plot function
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False,
            save="_furtherQC.pdf",
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="pct_counts_mt",
            show=False,
            save="_furtherQC_mt.png",
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="n_genes_by_counts",
            show=False,
            save="_furtherQC_gene.png",
        )

        print(f"Perform the further filtering process. The criteria:")
        print(
            f"% mitochondria: {pct_counts_mt_down_cutoff} - {pct_counts_mt_up_cutoff} %."
        )
        print(f"Total counts:  {total_counts_down_cutoff} - {total_counts_up_cutoff}.")
        print(
            f"Number of genes: {n_genes_by_count_down_cutoff} - {n_genes_by_count_up_cutoff} genes."
        )
        print(f"Cells beyond those criteria are excluded.")
        print(
            f"After the further filtering process, your data set contains {adata.n_obs} cells and {adata.n_vars} genes."
        )

        return adata


def make_analysis(adata):
    """generate analysis and plot"""
    pass


def save_csv():
    """ToDo: write results into csv"""
    pass


## main
if __name__ == "__main__":

    paths = get_input()
    exp = SCAnalysis(paths)
    adata = exp.load_data(paths)
    adata_filtered = exp.filter_data(adata)
    make_analysis(adata_filtered)

    save_csv()
