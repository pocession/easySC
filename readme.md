# easySC.py

## Description

This script performs a basic analysis of single-cell RNA-seq data using the scanpy package. It only supports 10X format at this moment.

## Dependencies

The following packages are required to run this code:

- `argparse`
- `pathlib`
- `numpy`
- `pandas`
- `scanpy`

You can run the following commands to run this script:

`pip install argparse pathlib numpy pandas scanpy`

Or install all packages through `requirements.txt`:

`pip3 install -r requirements.txt`

## Downloading codes

You can download the code either by cloning the repository from GitHub or by downloading the code as a ZIP file.

## Usage

- `usage: python3 src/easySC.py [-h] --data_dir DATA [--min_genes MIN_GENES] [--min_cells MIN_CELLS] [--max_cutoff MAX_CUTOFF] [--min_cutoff MIN_CUTOFF]`
- `--data_dir` `-d`: a folder contains three files: `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`.
- `--min_genes`: specify the minimum genes per cell (default = 200).
- `--min_cells`: keep genes that are expressed in the minimum number of cells (default = 3).
- `--max_cutoff`: keep cells with parameters under this percentile (default = 0.75).
- `--min_cutoff`: keep cells with parameters above this percentile (default = 0.25).
- `--fig_format`: set the file format of output figures (default = .png).

## Examples

- `python src/easySC.py --data_dir /path/to/data --min_genes 250 --min_cells 6`

This will perform a basic analysis of your single-cell RNA-seq data, keeping genes that are expressed in at least 6 cells and cells that have at least 250 expressed genes. Also keep cells with total count and gene count between the 25th and 75th percentile. The results will be saved to a a figures directory in png format as well as in results directory in csv format.

## Updates

- 2023-3-12: The first version is ready for the world.

## Discussion

- Please use `issue` in github to create issue and discussion.

## Test data

- `wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_raw_feature_bc_matrix.tar.gz`

## Unit tests

Use pytest for test.
`python3 -m pytest`

## Workflow

- ![workflow](easySC_workflow.png)
