# easySC.py

## Description

This app aims for producing preliminary profiles of single cell data. It only supports 10X format at this moment.

## Usage

- `usage: python3 src/easySC.py [-h] --data_dir DATA [--min_genes MIN_GENES] [--min_cells MIN_CELLS] [--max_cutoff MAX_CUTOFF] [--min_cutoff MIN_CUTOFF]`
- `--data_dir` `-d`: a folder contains three files: `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`.
- `--min_genes`: specify the minimum genes per cell (default = 200).
- `--min_cells`: keep genes that are expressed in the minimum number of cells (default = 3).
- `--max_cutoff`: keep cells with parameters under this percentile (default = 0.75).
- `--min_cutoff`: keep cells with parameters above this percentile (default = 0.25).
- `--fig_format`: set the file format of output figures (default = .png).

## Updates

## Discussion

- Please use `issue` in github to create issue and discussion.

## Test data

- `wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_raw_feature_bc_matrix.tar.gz`

## Unit tests

Use pytest for test.
`python3 -m pytest`

## Workflow

- ![workflow](easySC_workflow.png)
