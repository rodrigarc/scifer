#' run IgDiscover for IgBlast
#'
#' @param folder_sequences Folder containing all the sanger sequencing abi/ab1 files on subfolders. Each subfolder should have have a identifiable name, matching name with fcs data. eg. "E18_01", "E23_06". The first characters of the ab1 file name should be the well location. eg. "A1-sequence1.ab1", "F8_sequence-igg.ab1"
#'
#' @return List containing two items:
#' * summaries: contains all the summary results from the processed abi files,
#' * quality_scores: contains all the Phred quality score for each position.
#' @importFrom parallel detectCores mclapply
#' @importFrom sangerseqR read.abif
#' @importFrom DECIPHER ConsensusSequence
#' @rawNamespace import(Biostrings, except=c(collapse, union, intersect, setdiff, setequal))
#' @examples
#' sf <- summarise_quality(
#'     folder_sequences = system.file("extdata/sorted_sangerseq",
#'                                     package = "scifer"),
#'     secondary.peak.ratio = 0.33,
#'     trim.cutoff = 0.01,
#'     processor = 1
#' )
#' @export
library(scifer)
library(dplyr)

 sf <- summarise_quality(folder_sequences = system.file("extdata/sorted_sangerseq",
                                 package = "scifer"),
                         secondary.peak.ratio = 0.33,
                         trim.cutoff = 0.01,
                         processor = 4)
 pathnames <- as.character(sf$summaries$file.path)
 sangerseqlisted <- sapply(pathnames, sangerseqR::readsangerseq)
 sangerbasecall.string <- sapply(sangerseqlisted, sangerseqR::primarySeq,
                                 string=TRUE
 )
 sangerbasecall.string <- sangerbasecall.string %>%
   data.frame() %>%
   tibble::rownames_to_column()
 names(sangerbasecall.string)[names(sangerbasecall.string) == "rowname"] <- "file.path"
 names(sangerbasecall.string)[names(sangerbasecall.string) == "."] <- "sequence"
sangerbasecall.string$sequence

reticulate::conda_list()

library(reticulate)
use_condaenv("igxplore")
source_python("R/igblastwrap.py")
py_run_string("igblastwrap.py")



