# RepertoiR: Single-cell Immunglobulin Filtering of Sanger sequences
[![R build
status](https://github.com/rodrigarc/RepertoiR/workflows/R-CMD-check/badge.svg)](https://github.com/rodrigarc/RepertoiR/actions)

### Integrating index single-cell sorted files with Sanger sequencing per plates



Version: 0.1.0  

Author and Maintainer: Rodrigo Arcoverde  

Description: Have you ever index sorted cells in a 96 or 384-well plate and then sequenced using Sanger sequencing? If so, you probably had some strungle to either check the chromatograms of each cell sequenced manually, or identify which cell was sorted where after sequencing the plate. Scifer was developed to solve this issue by performing basic quality control of Sanger sequences and merging flow cytometry data from probed single-cell sorted B cells with sequencing data. RepertoiR can export summary tables, fasta files, chromatograms for visual inspection, and generate an interactive report.

## Installation

Before the installation you should first install the following packages from Bioconductor:

```r
# just copy and run the following to check for the needed packages and install them if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DECIPHER","sangerseqR","ape"))
```

To install RepertoiR directly from GitHub, run this:

```r
# just copy and run the following to install RepetoiR

if (!require("devtools"))
install.packages("devtools")
devtools::install_github("rodrigarc/RepertoiR")
```

## Generating reports and output files

The ```quality_report``` report function will look recursively on the selected folder for ab1 files for quality filtering. Check the ```test_dataset``` folder to be sure how to organize your sequences. Your ab1 files should start with the well number, and each plate should be under a different folder named accordingly to ID and plate. For example, ```E18_C1``` where E18 is the ```ID```, ```C1``` is the plate. Finally, repeated sequencing for the same sequences can also be filtered and the best sequence will be selected. For that, just add the repeated sequences in a folder with the same name as the original folder but add ```_R``` at the end, for example ```E18_C1_R```. 

The function should be run as the following:

```r
# You can use the following arguments 
# `data_folder` Full file directory for searching all ab1 files in a recursive search method (all files including files in subfolders)
# `outputfile` Output file name for the report generation
# `output_dir` Output directory for all the different output files that are generated during the report
# `processors` Number of processors to use, or NULL (the default) for all available processors
# `plot_chromatogram` Logical argument, TRUE or FALSE, to indicate if chromatograms should be plotted or not. Default is FALSE.
# `folder_path_fcs`  Full file directory for searching all flow cytometry index files (.fcs) in a recursive search method (all files including files in subfolders)


  
quality_report(data_folder = "full_path_to/test/sanger_sequences_folder",
               outputfile = "HC-QC-report.html",
               output_dir = "full_path_to/test_output_folder",
               processors = 1,
               plot_chromatogram = FALSE,
               folder_path_fcs = "full_path_to/test/fcs_index_files_folder)
               

```

This function will generate: 
 - general quality report (.html)
 - individualized report per ID (.html)
 - table containing the good quality sequences filtered in and their QC (.csv)
 - FASTA file containing all the sequences
 - FASTA files per ID 
 - Chromatogram of the approximmate CDR3 region if secondary peaks are detected (function not yet updated)  


## Special cases:

If you want to check the scores and filtering in a ```R object``` or do it for a single sequence, the package also contains functions for that which are called during ```quality_report``` function. To get a summary in a dataframe format of the quality of every sequence in the folder, you can use the following function. 

```r
#for many abi_files 
abi_files <- summarise_quality("/path/to/folder/with/ab1/files", processors = 1)
View(abi_files$summaries)

#for one single abi_files
abi_file <- summarise_abi_file("/path/to/ab1/file", processors = 1)
```

### FAQ
* I have an Mac with a M1 processer, how can I install the packages?
        
    First, install xcode:
        ```
        su - user - myyser -c 'xcode-select --install'
        ```

    Then install devtools, decipher, sangerseq, and ape.

        ```r
        install. packages("devtools")

        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

        BiocManager::install(c("DECIPHER","sangerseqR","ape"), force = TRUE)
        ```
