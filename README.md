# RepertoiR

### B cell receptor repertoire package development  



Version: 0.1.0  

Author and Maintainer: Rodrigo Arcoverde  

Description: RepertoiR is being developed for quality control of Sanger sequences, merging single-cell sorted specificity with sequencing data, and analysing/plotting repertoire data for exploratory analysis.  

## Installation

Before the installation you should first install the following packages from Bioconductor:

```r
# just copy and run the following to check for the needed packages and install them if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DECIPHER","sangerseqR","ape"))
```

To install directly from GitHub, run this:

```r
# just copy and run the following to install RepetoiR

if (!require("devtools"))
install.packages("devtools")
devtools::install_github("rodrigarc/RepertoiR")
```

## Generating reports and output files

The ```quality.report``` report function will look recursively on the selected folder for ab1 files for quality filtering. Check the ```test_dataset``` folder to be sure how to organize your sequences. Your ab1 files should start with the well number, and each plate should be under a different folder named accordingly to ID, animal, plate, chain. For example, ```E18_C1_HC``` where E18 is the ```ID```, ```C1``` is the plate, and ```HC``` stands for heavy chain. Finally, repeated sequencing for the same sequences can also be filtered and the best sequence will be selected. For that, just add the repeated sequences in a folder with the same name as the original folder but add ```_R``` at the end, for example ```E18_C1_HC_R```. 

The function should be run as the following:

```r
# data_folder is the root folder containing the ab1 file, be sure to add a valid location. 
## You can check that by running summarise.quality(path/to/ab1/folder) and seeing if you can detect ab1 files
# outputfile is the name for the report
# output_dir is the directory to add the files
# processors are the number of cores available in your computer for multithreading, 
#leave as 1 to run as single-thread, or NULL to automatically detect the number of cores
  
quality.report(data_folder = "full_path_to/test/teste_dataset/",
               outputfile = "HC-QC-report.html",
               output_dir = "full_path_to/test_output_folder",
               processors = 1)

```

This function will generate: 
 - general quality report (.html)
 - individualized report per ID (.html)
 - table containing the good quality sequences filtered in and their QC (.csv)
 - FASTA file containing all the sequences
 - FASTA files per ID 
 - Chromatogram of the approximmate CDR3 region if secondary peaks are detected (function not yet updated)  


## Special cases:

If you want to check the scores and filtering in a ```R object``` or do it for a single sequence, the package also contains functions for that which are called during ```quality.report``` function. To get a summary in a dataframe format of the quality of every sequence in the folder, you can use the following function. 

```r
#for many abi_files 
abi_files <- summarise.quality("/path/to/folder/with/ab1/files", processors = 1)
View(abi_files$summaries)

#for one single abi_files
abi_file <- summarise.abi.file("/path/to/ab1/file", processors = 1)
```

