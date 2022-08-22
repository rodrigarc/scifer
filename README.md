# scifer: Single-cell Immunoglobulin Filtering of Sanger sequences
<!-- badges: start -->
[![R build
status](https://github.com/rodrigarc/scifer/workflows/R-CMD-check/badge.svg)](https://github.com/rodrigarc/scifer/actions)
[![R-CMD-check-bioc](https://github.com/rodrigarc/scifer/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/rodrigarc/scifer/actions)
<!-- badges: end -->

### Integrating index single-cell sorted with Sanger sequencing raw files


Version: 0.99.0  

Author and Maintainer: Rodrigo Arcoverde Cerveira

Description: Have you ever index sorted cells in a 96 or 384-well plate and then sequenced using Sanger sequencing? If so, you probably had some struggle to either check the chromatograms of each cell sequenced manually, or identify which cell was sorted where after sequencing the plate. `scifer` was developed to solve this issue by performing basic quality control of Sanger sequences and merging flow cytometry data from probed single-cell sorted cells with sequencing data. `scifer` can export summary tables, fasta files, chromatograms for visual inspection, and generate a html report.
This package was developed focused in B cell receptor sequences from antigen-specific single-cell sorted B cells, however it is can be highly customizable to other type of single-cell sorted sequences.

## Installation

Before the installation you should first install the following packages from Bioconductor:

```r
# just copy and run the following to check for the needed packages and install them if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DECIPHER","sangerseqR","ape"))
```

To install scifer directly from GitHub, run this:

```r
# just copy and run the following to install scifer

if (!require("devtools"))
install.packages("devtools")
devtools::install_github("rodrigarc/scifer")
```

## Processing flow cytometry index data

If you want to merge the flow cytometry data with the sanger sequencing data, aka merging metadata from `.fcs` with `.ab1` files. Please first check the folder structure of your project. You should have 2 folders:
  1. Index flow cytometry data (with `.fcs` index files): this folder should contain all the index files that you want to process. It searchers recursively, so all the `.fcs ` within subfolders will also be used. The index file names should have the `ID` and the `plate` named, so it can be later on used for merging with the sanger data. 
  *eg. "~/Desktop/analysis/fcs_files", which contains all the files named such as "E18_02.fcs". where "E18" stands for the sample ID and "02" for the plate number.*
  
  2. Sanger sequencing ab1 file folders: this folder should contain all the sanger sequences with the `.ab1` extension. It also searchers recursively, so all the `.ab1` files within subfolders will also be used. In this folder you should have each sample ID divided into subfolders named with the sample ID and plate name. This is **CRUCIAL** since it will be used for matching the two different datasets. For the above example, you would have a subfolder named `E18_02` and within that folder all the `.ab1` files. The first characters separated by `_` or `-` should be the well plate position, which is usually how facilities name the results. 
   *eg. `A1_3_IgG_Inner.ab1`, where `A1` is the well position.*

Example of data folder structure:
 ```
project_folder
│
└─── results_output
│
└─── fcs_files
│         │ E18_01.fcs
│         │ E18_02.fcs
│         │ E24_06.fcs
│   
└───sanger_sequences
        └─── E18_01
        │       │ A1_3_IgG_Inner.ab1
        │       │ A2_3_IgG_Inner.ab1
        │       │ A3_3_IgG_Inner.ab1
        │       │ ...
        └─── E18_02
        │       │ A1_3_IgG_Inner.ab1
        │       │ H1_3_IgG_Inner.ab1
        │       │ ...
        └─── E24_06
                │ A1_3_IgG_Inner.ab1
                │ G12_IgG_Inner.ab1
                │...
```

After adjusting your data following this folder structure, you can first run `fcs_processing()` function to see your index sorted files merged and their fluorescence for each channel. It also plots a common flow plot to facilitate the threshold value for each probe you have used.

```r
library(scifer)

fcs_data <- fcs_processing(folder_path = "~/Desktop/project_folder/fcs_files",
                           compensation = TRUE, plate_wells = 96,
                           probe1 = "Pre.F", probe2 = "Post.F",
                           posvalue_probe1 = 600, posvalue_probe2 = 400))

```

## Inspecting sanger sequencing quality of single-cell sorts

After identifying which probes and threshold you wish to use for selecting what is positive and which channel name you should use. You can also inspect the sequencing data. This process takes longer than the previous steps depending on the number of sequences to be analyzed.

You can retrieve a `table` for inspection of all the sequences, and checking the names using the following function:

```r
summary_sanger_data <- summarise_quality(folder_sequences = "~/Desktop/project_folder/sanger_sequences",
                                         secondary.peak.ratio = 0.33,
                                         trim.cutoff = 0.01,
                                         processors = NULL # You can change the number of processors manually, 
                                         # or let the function decides how many are detected in your machine. More processors used, faster the results.
                                         )
```
You can check all the file names and quality measurements of each sequence, it can help you decide which quality thresholds you want to customize if needed.


## Generating reports and output files

The `quality_report()` function will use both of the above functions and others to generate different outputs. 
This function will generate: 
 - general quality report (.html)
 - individualized report per ID (.html)
 - table containing the good quality sequences filtered in and their QC (.csv)
 - `fasta` file containing all the sequences
 - `fasta` files per ID 
 - Chromatogram of the approximate CDR3 region if secondary peaks are detected 

The function should be run as the following:

```r
quality_report(folder_sequences = "~/Desktop/project_folder/sanger_sequences",
               outputfile = "QC_report.html",
               output_dir = "~/Desktop/project_folder/results_output",
               folder_path_fcs = "~/Desktop/project_folder/fcs_files",
               probe1 = "PE.Cy7.A", probe2 = "Alexa.Fluor.700.A",
               posvalue_probe1 = 600, posvalue_probe2 = 400)
               
# You can use the default arguments or completely customize the thresholds according to your dataset or desire, the following parameters can be customized:
# `plot_chromatogram` Logical argument, TRUE or FALSE, to indicate if chromatograms should be plotted or not. Default is FALSE
# `raw_length` Minimum sequence length for filtering. Default is 400 for B cell receptors
# `trim_start` Starting position where the sequence should start to have a good base call accuracy. Default is 50 for B cell receptors
# `trim_finish` Last position where the sequence should have a good base call accuracy. Default is 409 for B cell receptors
# `trimmed_mean_quality` Minimum Phred quality score expected for an average sequence. Default is 30, which means average of 99.9\% base call accuracy
# `compensation` Logical argument, TRUE or FALSE, to indicate if the index files were compensated or not. If TRUE, it will apply its compensation prior assigning specificities
# `plate_wells` Type of plate used for single-cell sorting. eg. "96" or "384", default is "96"
# `probe1` Name of the first channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
# probe2 Name of the second channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
# `posvalue_probe1` Threshold used for fluorescence intensities to be considered as positive for the first probe
# `posvalue_probe2` Threshold used for fluorescence intensities to be considered as positive for the second probe
# `cdr3_start` Expected CDR3 starting position, that depends on your primer set. Default is position 100
# `cdr3_end` Expected CDR3 end position, that depends on your primer set. Default is position 150

```


## Special cases:

#### Inspecting a single ab1 file

If you want to inspect only one ab1 file, maybe in case the `summarise_quality()` is giving your errors, you can do that by using the `summarise_abi_file()` function. 

eg.
```r
summarise_abi_file("~/Desktop/project_folder/sanger_sequencing/A1_3_IgG_Inner.ab1")
```

#### Saving fasta file

If you want to save manually a `data.frame` or a vector to a `.fasta file`, you can use:

```r
df <- data.frame(sequence_name = c("seq1", "seq2"), sequence = c("CGATGCAT", "ATCGAGCTG"))
df_to_fasta(sequence_name = df$sequence_name,
            sequence_strings = df$sequence,
            file_name = "my_fasta.fasta",
            output_dir = "~/Desktop/results_output")
```

## FAQ
* I have a Mac with a M1 processor and I am getting errors while installing this package. How can I install it?
        
First, install xcode:

```
    #open your teminal and type
    xcode-select --install
    
    #if the above does not work due to lack of user permissions, you can use use a specific admin user to call the function:
    su - user - my_user -c 'xcode-select --install'  
    
```

Then install devtools, decipher, sangerseq, and ape.
   

```r
    install. packages("devtools")

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install(c("DECIPHER","sangerseqR","ape"), force = TRUE)
```
