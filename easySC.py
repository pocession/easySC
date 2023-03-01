## library
import sys
import os.path

## functions

def examine_argv():
    ## examine the arguments
    if len(sys.argv) < 3:
        print("Please provide both input and output")
    
    # Iterate through each path and check if it exists
    paths = sys.argv[1:]
    for path in paths:
        if os.path.exists(path):
            print(f"The path '{path}' exists!")
        else:
            print(f"The path '{path}' does not exist.")

    if len(paths) >=2:
        inPath = paths[0]
        outPath = paths[1]

    # Check input files 
    ## Check if there is any files in input folder
    if os.listdir(inPath) == []:
        print("No files found in the input directory.")
    else:
        print("Some files found in the input directory.")

    barcodes = os.path.join(inPath,"barcodes.tsv")
    features = os.path.join(inPath,"features.tsv")
    matrix = os.path.join(inPath,"matrix.mtx")

    if not os.path.isfile(barcodes):
        print("Please provide the barcode file named barcodes.tsv!")
    
    if not os.path.isfile(features):
        print("Please provide feature (gene) file named features.tsv")

    if not os.path.isfile(matrix):
        print("Please provide the count matrix named matrix.mtx")

def get_input():
    print("Get input")

def filter_data():
    ## ToDo: perform basic filtering, umi/gene, gene/cell
    print("Filter data")

def make_analysis():
    ## generate analysis and plot
    print("make analysis")

def save_csv():
    ## ToDo: write results into csv
    print("Save CSV")

def save_plot():
    ## save plots
    with open(outFile,'w') as o:
        for line in processedLines:
            o.write(line)


## main
examine_argv()
""" input = get_input()
filtered = filter_data(input)
make_analysis(filtered)

save_csv()
save_plot() """