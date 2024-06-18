#' run IgDiscover for IgBlast

library(reticulate)

#check which environments are installed
conda_list()
# conda env created through the terminal
# conda create -n igblast-wrap python=3.9 dnaio igblast
# load conda environment within R
use_condaenv("igblast-wrap")
# run igblast
source_python("R/igblastwrap.py")
output <- system("python R/igblastwrap.py 'inst/extdata/test_fasta/KIMDB_rm' 'inst/extdata/test_fasta/test_igblast.txt'",
                 intern = TRUE)
# get airr formated table results
airr_results <- read.table(text = output, header = TRUE, sep = "\t")
airr_results






