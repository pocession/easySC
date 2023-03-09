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

## Functions
### Get input data

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

            print(f"missing file: {f}")


    return args

### Save to csv files
def save_csv(df,str):
    """ToDo: write results into csv"""
    compression = 'zip'
    filename = (str,'csv',compression)
    filename = ".".join(filename)
    filepath = Path('./results/',filename)
    compression_opts = dict(method='zip')
    filepath.parent.mkdir(parents=True, exist_ok=True) 
    df.to_csv(filepath,sep=",",compression = compression_opts)


class SCAnalysis:
    def __init__(self, args=list()):
        self.args = args

        self.adata = None

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

        self.adata = sc.read_10x_mtx(
            args.data_dir,  # the directory with the `.mtx` file
            var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
            cache=True,
        )
        self.adata.var_names_make_unique()  # necessary because we use gene symbols as var_names
        print(f"Your data set contains {self.adata.n_obs} cells and {self.adata.n_vars} genes.")


    def filter_data(self):
        """
        ToDo: perform filtering and produce QC plots
        """

        # p.data,p.min_genes,p.min_cells,p.up_cutoff,p.min_cutoff

        min_genes = args.min_genes
        min_cells = args.min_cells
        up_cutoff = args.max_cutoff
        down_cutoff = args.min_cutoff
        adata = self.adata

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


        print(
            f"After the further filtering process, your data set contains {adata.n_obs} cells and {adata.n_vars} genes."
        )

        self.adata = adata

    def plot_filtered_data(self):
        # Draw figures for raw data and save them.
        # The draw function comes from scanpy itself.
        # # The out plots will be save in "./figures" automatically.


        adata = self.adata

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

            save="_preliminaryQC_mt.pdf",

        )
        sc.pl.scatter(
            self.data,
            x="total_counts",
            y="n_genes_by_counts",
            show=False,

            save="_preliminaryQC_gene.pdf",

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

            save="_furtherQC_mt.pdf",

        )
        sc.pl.scatter(
            self.data,
            x="total_counts",
            y="n_genes_by_counts",
            show=False,

            save="_furtherQC_gene.pdf",
        )


    def make_analysis(self):
        """generate analysis and plot"""

        adata = self.adata

        # Total-count normalize
        sc.pp.normalize_total(adata, target_sum=1e4)

        # Logarithmize the data:
        sc.pp.log1p(adata)

        # Identify highly-variable genes.
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

        # freeze the data into raw attribute
        adata.raw = adata

        ## Save the high variant gene plot
        sc.pl.highly_variable_genes(adata,show=False,save="_highlyVariableGenes.pdf")

        ## Filter the data based on highly variable genes
        adata_hvg = self.adata[:,adata.var.highly_variable]

        ## Not sure what the below two steps are doing
        ### Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
        ### Scale each gene to unit variance. Clip values exceeding standard deviation 10.

        # sc.pp.regress_out(adata_hvg, ['total_counts', 'pct_counts_mt'])
        # sc.pp.scale(adata_hvg, max_value=10)

        df_hvg = pd.DataFrame.sparse.from_spmatrix(adata_hvg.X,index=adata_hvg.obs_names,columns=adata_hvg.var_names)
        save_csv(df_hvg,"hvg")


        #umap
        sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata_hvg)
        sc.tl.leiden(adata_hvg)
        sc.pl.umap(adata_hvg,color=['leiden'],use_raw=False,show=False,save="_leiden_pdf")


## main
if __name__ == "__main__":

    args = get_input()
    exp = SCAnalysis(args)
    exp.load_data()
    exp.filter_data()
    exp.plot_filtered_data()
    exp.make_analysis()
