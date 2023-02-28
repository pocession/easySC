## library

## functions

def get_input():
    ## ToDO: read data and check input

def filter_data():
    ## ToDo: perform basic filtering, umi/gene, gene/cell

def make_analysis():
    ## generate analysis and plot

def save_csv():
    ## ToDo: write results into csv

def save_plot():
    ## save plots


## main
input = get_input()
filtered = filter_data(input)
make_analysis(filtered)

save_csv()
save_plot()