# scifer: Single-cell Immunoglobulin Filtering of Sanger sequences <img src="man/figures/logo.png" align="right" width="180" />
<!-- badges: start -->
[![R build
status](https://github.com/rodrigarc/scifer/workflows/master_r-cmd-check/badge.svg)](https://github.com/rodrigarc/scifer/actions)
[![codecov](https://codecov.io/gh/rodrigarc/scifer/branch/r-cmd-check/graph/badge.svg)](https://codecov.io/gh/rodrigarc/scifer)
<!-- badges: end -->

### Integrating index single-cell sorted with Sanger sequencing raw files


Version: 1.8.0  

Author and Maintainer: Rodrigo Arcoverde Cerveira

Description: Have you ever index sorted cells in a 96 or 384-well plate and then sequenced using Sanger sequencing? If so, you probably had some struggle to either check the electropherograms of each cell sequenced manually, or identify which cell was sorted where after sequencing the plate. `scifer` was developed to solve this issue by performing basic quality control of Sanger sequences and merging flow cytometry data from probed single-cell sorted cells with sequencing data. `scifer` can export summary tables, fasta files, electropherograms for visual inspection, and generate a html report.
This package was developed focused in B cell receptor sequences from antigen-specific single-cell sorted B cells, however it is highly customizable to other type of single-cell sorted sanger sequences. 

* `scifer` has been primarily tested with R versions >= 4.2.3 and Bioconductor versions >= 3.16. It is recommended to use these or newer versions of R and Bioconductor for best results. 

## Installation

Before the installation, you should first install the following packages from Bioconductor:

```r
# just copy and run the following to check for the needed packages and install them if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DECIPHER","sangerseqR","ape", "BiocStyle"))
```

If versions of `DECIPHER`, `sangerseqR`, `ape`, and `BiocStyle` are already installed, you may need to set `force` = `TRUE` to update them.

```r
# just copy and run the following to install the needed packages

BiocManager::install(c("DECIPHER","sangerseqR","ape", "BiocStyle"), force = TRUE)

```
`scifer` with `ìgblast` functionality is currently available for MacOS-64 and windows platforms. 

### Windows

In order to enable `ìgblast` on windows the user must first download and run `ncbi-igblast-1.22.0-win64.exe` from the NCBI igblast website [here](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/). Do not change the default installation location. If R-studio is open, please restart R-studio for the changes to take effect. 

### MacOS-arm64 (M1/M2/M3 - Apple silicon processors)

`igblast` functionality is not currently available for MacOS-arm64 platforms through `bioconda`. However, `scifer` can still use the function `igblast()` on MacOS-arm64 platforms by using the following steps: 

1. Download version R-4.2.3 from CRAN: [here](https://cran.r-project.org/bin/macosx/base/R-4.2.3.pkg). Installing this version of R will allow `scifer` to use the `igblast()` function on MacOS-arm64 platforms once installed and should automatically become the default version of R when opening R-studio.

2. Download and install the following packages before installing `scifer`: 
```r
install.packages(c("dplyr", "rmarkdown", "data.table", "plyr", "knitr", "stringr", "ggplot2", "gridExtra", "tibble", "scales", "rlang", "reticulate", "here"))

BiocManager::install(c("basilisk", "basilisk.utils", "flowCore", "BiocStrings"))
```
If asked to compile from source, please select `yes`. 

### Bioconductor or GitHub scifer versions

When installing `scifer` you should choose between installation from GitHub or from Bioconductor. The GitHub version is the most recent one with updated developmental features, changes, and bug fixes. The GitHub version is used for testing new features and bug fixes before they are submitted to Bioconductor. The Bioconductor version is the most stable one and compliant with Bioconductor's submission process. This version is updated every 6 months.

To install scifer directly from GitHub, run this:

```r
# just copy and run the following to install scifer

if (!require("devtools"))
install.packages("devtools")
devtools::install_github("rodrigarc/scifer")
```

OR,

To install scifer directly from Bioconductor, run this:

```r
# just copy and run the following to install scifer

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scifer")
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

# Select the folder location of the fcs files (flow cytometry index data), this example uses the data provided with `scifer`
directory_flowdata <- system.file("extdata/fcs_index_sorting",
                                   package = "scifer")
                                   
fcs_data <- fcs_processing(folder_path = directory_flowdata,
                           compensation = TRUE, plate_wells = 96,
                           probe1 = "Pre.F", probe2 = "Post.F",
                           posvalue_probe1 = 600, posvalue_probe2 = 400)

```

## Inspecting sanger sequencing quality of single-cell sorts

After identifying which probes and threshold you wish to use for selecting what is positive and which channel name you should use. You can also inspect the sequencing data. This process takes longer than the previous steps depending on the number of sequences to be analyzed.

You can retrieve a `table` for inspection of all the sequences, and checking the names using the following function:

```r
# Select the folder location of the fcs files (flow cytometry index data), this example uses the data provided with `scifer`
directory_sequences <- system.file("extdata/sorted_sangerseq",
                                    package = "scifer")
                                    
summary_sanger_data <- summarise_quality(folder_sequences = directory_sequences,
                                         secondary.peak.ratio = 0.33,
                                         trim.cutoff = 0.01,
                                         processor = NULL # This will try to use all your processors, but you can also change it to a desired number
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
# Select the folder location of the ab1 files (sanger sequences), this example uses the data provided with `scifer`
directory_sequences <- system.file("extdata/sorted_sangerseq",
                                    package = "scifer") 
# Select the folder location of the fcs files (flow cytometry index data), this example uses the data provided with `scifer`
directory_flowdata <- system.file("extdata/fcs_index_sorting",
                                   package = "scifer")
                                   
quality_report(folder_sequences = directory_sequences,
               outputfile = "QC_report.html",
               output_dir = "~/Desktop/project_folder/results_output",
               folder_path_fcs = directory_flowdata,
               probe1 = "PE.Cy7.A", probe2 = "Alexa.Fluor.700.A",
               posvalue_probe1 = 600, posvalue_probe2 = 400)
               
# You can use the default arguments or completely customize the thresholds according to your dataset or desire, the following parameters can be customized:
# `plot_chromatogram` Logical argument, TRUE or FALSE, to indicate if chromatograms should be plotted or not. Default is FALSE
# `raw_length` Minimum sequence length for filtering. Default is 400 for B cell receptors
# `trim_start` Starting position where the sequence should start to have a good base call accuracy. Default is 50 for B cell receptors
# `trim_finish` Last position where the sequence should have a good base call accuracy. Default is 409 for B cell receptors
# `trimmed_mean_quality` Minimum Phred quality score expected for an average sequence. Default is 30, which means average of 99.9\% base call accuracy
# `compensation` Logical argument, TRUE or FALSE, to indicate if the index files were compensated or not. If TRUE, it will apply its compensation prior to assigning specificities
# `plate_wells` Type of plate used for single-cell sorting. eg. "96" or "384", default is "96"
# `probe1` Name of the first channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
# probe2 Name of the second channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
# `posvalue_probe1` Threshold used for fluorescence intensities to be considered as positive for the first probe
# `posvalue_probe2` Threshold used for fluorescence intensities to be considered as positive for the second probe
# `cdr3_start` Expected CDR3 starting position, that depends on your primer set. Default is position 100
# `cdr3_end` Expected CDR3 end position, that depends on your primer set. Default is position 150

```
## Default parameters:

The default parameters for the quality control of sequences derived from a data set analyzing B cell receptor sequences from rhesus macaques (_Macaca_ _mulatta_) generated by [Ols et al., _Immunity_ 2023](https://doi.org/10.1016/j.immuni.2023.08.011). To generate this data set, primers from [Sundling et al., _Sci._ _Transl._ _Med_ 2013](https://doi.org/10.1126/scitranslmed.3003752) and [_J._ _Immunol._ _Methods_ 2013](https://doi.org/10.1016/j.jim.2012.09.003) were used. Those primers were developed to capture the entire V region from heavy and light chain from macaque B cell receptors, the mean length of those sequences was 460 basepairs (bp). The default parameters were also tested for a human B cell receptor dataset generated using primers from [Doria-Rose et al., _J._ _Virol._ 2016](https://doi.org/10.1128/jvi.01791-15). Additionally, we have tested scifer for gamma delta T cell receptor sequences from human and mice (_Mus_ _musculus_), for those it was required to change the default parameters due to its shorter amplicon length of 250 bp based on the primers used. The most important information is the length of your expected amplicon, with a different amplicon length you should change the default parameters, most importantly the `raw_length` and `trim_finish` parameters since they rely on your amplicon length. The default minimum mean Phred quality score equal to 30 is standard to achieve mean of 99.9% base call accuracy. Below you can find the defaults that are used by scifer and what you should change for other species or amplicons. These defaults can be used in other experiments or can be changed dependent on the users' needs. Below are defaults for B cell receptors for macaques/humans and T cell receptors for mice. Please see the vignette for more information on using `scifer`: [scifer walkthrough vignette](https://bioconductor.org/packages/release/bioc/vignettes/scifer/inst/doc/flow_cytometry_data.html)     

Default parameters for B cell receptors in macaques and humans as mentioned above: 

```r
# `raw_length` = 343 # based on your expected amplicon length
# `trim_start` = 50 # start of expected good quality base call position used for filtering sequences, based on your primer set
# `trim_finish` = 409 #  end of expected good quality base call position used for filtering sequences, based on your primer set
# `trimmed_mean_quality` = 30 # 99.9% base call accuracy
# `cdr3_start` = 100 # expected start of CDR3 position, based on your primer set
# `cdr3_end` = 150 # expected end of CDR3 position, based on your primer set

```

Default parameters for gamma delta T cell receptors in mice (_Mus musculus_): 

```r
# `raw_length` = 200 # based on your expected amplicon length
# `trim_finish` = 250 # used for filtering sequences, based on your primer set

```


Default parameters for when uncertain:

If you are uncertain of your amplicon length, first check the quality metrics of your sequences using `summarise_quality`. Save the results in an object to evaluate the quality metrics columns, especially `raw.length` and `trim.finish` columns, such as the mean, min, and max values. This will help you to identify which parameters you should change.

```r
quality_control_summary <- summarise_quality(folder_sequences = directory_sequences)

```

If you are still uncertain you can use a broad default parameters that will filter mostly on mean base call accuracy. To do that, change the `raw_length` to 0 and the `trim_finish` to 50. This will null the filter based on length but the base mean accuracy  of 99.9% will be maintained. You can also put `cdr3_start` and `cdr3_end` to 0 if there is no position that you are particularly interested in having with higher accuracy.

```r
quality_report(folder_sequences = directory_sequences,
              outputfile = "QC_report.html",
              output_dir = "~/Desktop/project_folder/results_output",
              raw_length = 0, trim_finish = 50,
              cdr3_start = 0, cdr3_end = 0)
```

Default for flow cytometry data:

Flow cytometry data should be examined to determine the optimal parameters based on your dataset and experiment. It is important to check your flow data to see how the data is being processed, if it is already compensated, and if the cells were probed. These conditions affect which thresholds you should use.

In order to see the available channels and thresholds use the function `fcs_processing()` to create an object `fcs_data`, after which you can see the available channels using: 

```{r channels}
   colnames(fcs_data)

```

The `specificity` column is based on your selected channels and their thresholds and uses these thresholds to name your sorted cells. The value used as input to `posvalue_probe1` and `posvalue_probe2` is the Mean Fluorescence Intensity (MFI) for a given channel. Probed samples would have a `Pre.F` single-positive, `Post.F` single-positive, double-positive cells named as`DP`, and double-negative cells named as `DN`. Which can be seen using: 

```{r specificity, message=FALSE, warning=FALSE}
   unique(fcs_data$specificity)

```

If cells were not probed for a specific antigen, you can just use the following thresholds and ignore the `specificity` column. This approach will add all of your cells, regardless of the detected fluorescence in a channel.

```r
 fcs_data <- fcs_processing( 
                 posvalue_probe1 = 0, 
                 posvalue_probe2 = 0
             )
```

The default for `compensation` is `FALSE` though often the samples were already compensated before sorting. Setting `compensation` to `TRUE` makes it so that compensation matrix within the index files will be already automatically applied.

```{r check_fcs_probe_compensated, message=FALSE, warning=FALSE}
fcs_data <- fcs_processing(
  folder_path = system.file("/extdata/fcs_index_sorting",package = "scifer"),
  compensation = TRUE, plate_wells = 96,
  probe1 = "Pre.F", probe2 = "Post.F",
  posvalue_probe1 = 600, posvalue_probe2 = 400)

fcs_plot(fcs_data)
```

At the moment scifer only accepts merging and integrating data with two channels (eg. 2 probes, one per channel). Furthermore, scifer does not do batch correction for flow cytometry data, thus different experiments should have their Median Fluorescence Intensity (MFI) normalized before running fcs_processing.

For more information please consult the vignette: [scifer walkthrough vignette](https://bioconductor.org/packages/release/bioc/vignettes/scifer/inst/doc/flow_cytometry_data.html)

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

* I have a Mac with a M1 processor and I am getting errors while installing this package. How can I fix it?
        
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

* I have a Windows computer and I am having problems to build `scifer`. What can I do?

You make sure that you have installed the most recent Rtools, which can be found [here](https://cran.r-project.org/bin/windows/Rtools/). After installing it, restart your R session and then try to build it again. 

* I have a Mac with a M1/M2/M3 processor and I have encountered an error running `igblast()`. What can I do?

An unexpected error may occur the first time running `ìgblast()` which fails to adequately set the environment variable for `python=3.9.19`. Running `igblast()` a second time creates the correct environment variable and the function runs as expected.

* I have a Mac with a M1/M2/M3 processor and am having issues with the following error: `mach-o file, but is an incompatible architecture (have 'arm64', need 'x86_64'))` when installing the package dependencies. What can I do?

This error is due to the fact that the package dependencies are not yet available for the arm64 architecture. To solve this issue, you must insure that the correct version of R is being used. This can be achieved by using RSwitch for updating the `R.framework/Versions/Current` directory following the instructions found here: [here](https://support.posit.co/hc/en-us/articles/200486138-Changing-R-versions-for-the-RStudio-Desktop-IDE). 

* When I run `igblast()` it takes too long. What can I do?

It is normal that the first `igblast()` call is slow. This function is a wrapper that relies on conda environments that are built during the first time you use it. After the environments are created, reusing the function will be faster since the environments are already built.

## Code of Conduct

Please note that the `scifer` project is released with a [Contributor Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

## Development tools

* Continuous code testing is possible thanks to [GitHub actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)  through `r BiocStyle::CRANpkg('usethis')`, `r BiocStyle::CRANpkg('remotes')`, and `r BiocStyle::CRANpkg('rcmdcheck')` customized to use [Bioconductor's docker containers](https://www.bioconductor.org/help/docker/) and `r BiocStyle::Biocpkg('BiocCheck')`.
* Code coverage assessment is possible thanks to [codecov](https://codecov.io/gh) and `r BiocStyle::CRANpkg('covr')`.
* The code is styled automatically thanks to `r BiocStyle::CRANpkg('styler')`.
* The documentation is formatted thanks to `r BiocStyle::CRANpkg('devtools')` and `r BiocStyle::CRANpkg('roxygen2')`.
