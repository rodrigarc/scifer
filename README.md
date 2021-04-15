# RepertoiR

### B cell receptor repertoire package development 

*Version: 0.1.0
*Author and Maintainer: Rodrigo Arcoverde

*Description: RepertoiR is being developed for quality control of sanger sequences, merging single-cell sorted specificity with sequencing data, and analysing/plotting repertoire data for exploratory analysis.  


## Installation

To install directly from GitHub, run this:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("rodrigarc/RepertoiR")
```

If the devtools-based approach does not work, you can download one of the built tar-balls from the builds directory and manually install the package from source by executing the following lines (replace DOWNLOAD_FOLDER with the absolute path to the tar-ball and VERSION with the version number of the downloaded package):

```r
install.packages("/DOWNLOAD_FOLDER/ImSpectR_VERSION.tar.gz",
                     repos = NULL,
                     type  = "source")
```
## How to use this package for processing single-cell sanger sequences data:

To get a summary in a dataframe format of the quality of every sequence in the folder, you can use the following function. 
```r
#for many abi_files 
abi_files <- summarise.quality("/path/to/folder/with/ab1/files", processors = 1)
View(abi_files$summaries)

#for one single abi_files
abi_file <- summarise.abi.file("/path/to/ab1/file", processors = 1)
```


