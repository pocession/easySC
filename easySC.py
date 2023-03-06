## library

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

print("..  loading packages ..")
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
print("...  loading packages complete!")

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
        "--data_dir", "-d", required=True, type=Path, help="specify the input path"
    )
    parser.add_argument(
        "--min_genes",
        required=False,
        type=int,
        default=200,
        help="specify the minimum genes per cell (default = 200)",
    )
    parser.add_argument(
        "--min_cells",
        required=False,
        type=int,
        default=3,
        help="keep genes that are expressed in the minimum number of cells (default = 3)",
    )
    parser.add_argument(
        "--max_cutoff",
        required=False,
        type=int,
        default=0.75,
        help="keep cells with parameters under this percentile (default = 0.75)",
    )
    parser.add_argument(
        "--min_cutoff",
        required=False,
        type=int,
        default=0.25,
        help="keep cells with parameters above this percentile (default = 0.25)",
    )

    args = parser.parse_args()

    # Check required input data files
    files = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]
    for f in files:
        f = args.data_dir / f
        if f.is_file():
            print(f"file found: {f}")
        else:
            print(f"File not found: {f}")

    return args


class SCAnalysis:
    def __init__(self, args=list()):
        self.args = args
        self.data = None

    def load_data(self):
        """
        Load the 3 data files
        Save to a data class object
        scanpy read all three files and put it in an adata object. Not sure if we still need a class object here.
        adata: https://anndata.readthedocs.io/en/latest/
        scanpy: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
        """

        # print(paths)
        print("..  reading input data..")
        self.data = sc.read_10x_mtx(
            self.args.data_dir,  # the directory with the `.mtx` file
            var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
            cache=True,
        )
        self.data.var_names_make_unique()  # necessary because we use gene symbols as var_names
        print(f"Your data set contains {self.data.n_obs} cells and {self.data.n_vars} genes.")

    def filter_data(self):
        """
        ToDo: perform filtering and produce QC plots
        """

        # p.data,p.min_genes,p.min_cells,p.up_cutoff,p.min_cutoff
        min_genes = self.args.min_genes
        min_cells = self.args.min_cells
        up_cutoff = self.args.max_cutoff
        down_cutoff = self.args.min_cutoff
        adata = self.data

        # Preliminary filtering process
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        # Generate QC metrics
        adata.var["mt"] = adata.var_names.str.startswith(
            "MT-"
        )  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        # Transform the objects of adata into a data frame
        df_obs = pd.DataFrame(adata.obs)

        # Determine the cutoff limits
        pct_counts_mt_up_cutoff = df_obs["pct_counts_mt"].quantile(up_cutoff)
        total_counts_up_cutoff = df_obs["total_counts"].quantile(up_cutoff)
        n_genes_by_count_up_cutoff = df_obs["n_genes_by_counts"].quantile(up_cutoff)

        # It is not necessary to filter cells based on the minimal % mitochondrial RNA
        # pct_counts_mt_down_cutoff = df_obs["pct_counts_mt"].quantile(down_cutoff)
        total_counts_down_cutoff = df_obs["total_counts"].quantile(down_cutoff)
        n_genes_by_count_down_cutoff = df_obs["n_genes_by_counts"].quantile(down_cutoff)

        # Further filtering process
        adata = adata[adata.obs.total_counts < total_counts_up_cutoff, :]
        adata = adata[adata.obs.n_genes_by_counts < n_genes_by_count_up_cutoff, :]
        adata = adata[adata.obs.pct_counts_mt < pct_counts_mt_up_cutoff, :]

        adata = adata[adata.obs.total_counts > total_counts_down_cutoff, :]
        adata = adata[adata.obs.n_genes_by_counts > n_genes_by_count_down_cutoff, :]
        # adata = adata[adata.obs.pct_counts_mt > pct_counts_mt_down_cutoff, :]
        print(f"The pipeline process data with the following criteria:")
        print(f"\tCells should contain more than {min_genes} expressed genes.")
        print(f"\tGenes should be expressed in more than {min_cells} cells.")
        print(
            f"\tCells should contain less than {pct_counts_mt_up_cutoff}% mitochondrial RNA (within {up_cutoff} percentile)."
        )
        print(
            f"\tCells should contain total counts with the following range:  {total_counts_down_cutoff} - {total_counts_up_cutoff} ({down_cutoff} - {up_cutoff} percentile)."
        )
        print(
            f"\tCells should contain genes within the following range: {n_genes_by_count_down_cutoff} - {n_genes_by_count_up_cutoff} (represents {down_cutoff} - {up_cutoff} percentile)."
        )

        self.data = adata

    def plot_filtered_data(self):
        # Draw figures for raw data and save them.
        # The draw function comes from scanpy itself.
        # # The out plots will be save in "./figures" automatically.

        ## This should be wrapped up into a plot function
        sc.pl.violin(
            self.data,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False,
            save="_preliminaryQC.pdf",
        )
        sc.pl.scatter(
            self.data,
            x="total_counts",
            y="pct_counts_mt",
            show=False,
            save="_preliminaryQC_mt.png",
        )
        sc.pl.scatter(
            self.data,
            x="total_counts",
            y="n_genes_by_counts",
            show=False,
            save="_preliminaryQC_gene.png",
        )

        ## This should be wrapped up into a plot function
        sc.pl.violin(
            self.data,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False,
            save="_furtherQC.pdf",
        )
        sc.pl.scatter(
            self.data,
            x="total_counts",
            y="pct_counts_mt",
            show=False,
            save="_furtherQC_mt.png",
        )
        sc.pl.scatter(
            self.data,
            x="total_counts",
            y="n_genes_by_counts",
            show=False,
            save="_furtherQC_gene.png",
        )

        print(
            f"Finished filtering. Data size: \n",
            "Cells: {self.data.n_obs} \n",
            "Genes: {self.data.n_vars} "
        )


    def make_analysis(self):
        """generate analysis and plot"""
        pass


    def save_csv(self):
        """ToDo: write results into csv"""
        pass


## main
if __name__ == "__main__":

    args = get_input()
    exp = SCAnalysis(args)
    exp.load_data()
    exp.filter_data()
    exp.plot_filtered_data()
    exp.make_analysis()
    exp.save_csv()
