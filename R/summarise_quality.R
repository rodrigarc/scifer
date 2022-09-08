#' Generate a summary table containing quality measurements from sanger sequencing abi files
#'
#' @param folder_sequences Folder containing all the sanger sequencing abi/ab1 files on subfolders. Each subfolder should have have a identifiable name, matching name with fcs data. eg. "E18_01", "E23_06". The first characters of the ab1 file name should be the well location. eg. "A1-sequence1.ab1", "F8_sequence-igg.ab1"
#' @param trim.cutoff Cutoff at which you consider a base to be bad. This works on a logarithmic scale, such that if you want to consider a score of 10 as bad, you set cutoff to 0.1; for 20 set it at 0.01; for 30 set it at 0.001; for 40 set it at 0.0001; and so on. Contiguous runs of bases below this quality will be removed from the start and end of the sequence. Given the high quality reads expected of most modern ABI sequencers, the defualt is 0.0001.
#' @param secondary.peak.ratio Ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated, while those below the ratio are not.
#' @param processors Number of processors to use, or NULL (the default) for all available processors
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
summarise_quality <-
    function(folder_sequences = "input_folder",
    trim.cutoff = 0.01,
    secondary.peak.ratio = 0.33,
    processors = NULL) {
        get_processors <- function(processors) {
            if (Sys.info()["sysname"] == "Windows") {
                ## mclapply is not supported on windows
                return(1)
            }
            if (is.null(processors)) {
                processors <- detectCores(all.tests = FALSE, logical = FALSE)
            }
            return(processors)
        }
        processors <- get_processors(processors)
        message("Looking for .ab1 files...")
        if (!dir.exists(folder_sequences)) {
            stop("Folder containing sequences does not exist.")
        } else {
            abi.fnames <- list.files(folder_sequences,
                pattern = "\\.ab1$",
                full.names = TRUE, recursive = TRUE
            )
            message(sprintf(("Found %d .ab1 files..."), length(abi.fnames)))
            message("Loading reads...")
            abi.seqs <- mclapply(abi.fnames,
                                 sangerseqR::read.abif,
                                 mc.cores = processors)
            message("Calculating read summaries...")
            ## Create a data.frame of summaries of all the files
            summaries.dat <- mclapply(abi.seqs,
                summarise_abi_file,
                trim.cutoff = trim.cutoff,
                secondary.peak.ratio = secondary.peak.ratio,
                mc.cores = processors
            )
            message("Cleaning up")
            summaries <- mclapply(summaries.dat,
                                  function(x) x[["summary"]],
                                  mc.cores = processors)
            summaries <- do.call(rbind, summaries)

            folder.names <- basename(dirname(abi.fnames))
            file.names <- basename(abi.fnames)

            summaries <- cbind.data.frame(
                "file.path" = as.character(abi.fnames),
                "folder.name" = as.character(folder.names),
                "file.name" = file.names, summaries, stringsAsFactors = FALSE
            )
            qual_scores <- mclapply(summaries.dat,
                                    function(x) x[["quality_score"]],
                                    mc.cores = processors)
            names(qual_scores) <- as.character(abi.fnames)

            return(list("summaries" = summaries,
                        "quality_scores" = qual_scores))
        }
    }
