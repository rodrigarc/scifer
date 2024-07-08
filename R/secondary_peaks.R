#' Check for secondary peaks in a sangerseq object
#'
#' This function finds and reports secondary peaks in a sangerseq object.
#' It returns a table of secondary peaks, and optionally saves an annotated chromatogram and a csv file of the peak locations.
#'
#' @param s a sangerseq s4 object from the sangerseqR package
#' @param ratio Ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not.
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param file.prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' @param processors Number of processors to use, or NULL (the default) for all available processors
#'
#' @return A list with two elements:
#'          \enumerate{
#'              \item {secondary.peaks}: a data frame with one row per secondary peak above the ratio, and three columns: "position" is the position of the secondary peak relative to the primary sequence; "primary.basecall" is the primary base call; "secondary.basecall" is the secondary basecall. \cr
#'              \item {read}: the input sangerseq s4 object after having the makeBaseCalls() function from sangerseqR applied to it. This re-calls the primary and secondary bases in the sequence, and resets a lot of the internal data.
#'          }
#'
#' @rawNamespace import(Biostrings, except=c(collapse, union, intersect, setdiff, setequal))
#' @importFrom sangerseqR primarySeq makeBaseCalls secondarySeq chromatogram
#' @importFrom DECIPHER AlignSeqs
#' @importFrom stringr str_locate_all
#' @importFrom methods is
#'
#' @examples
#' ## Read abif using sangerseqR package
#' s4_sangerseq <- sangerseqR::readsangerseq(
#'     system.file("/extdata/sorted_sangerseq/E18_C1/A1_3_IgG_Inner.ab1",
#'         package = "scifer"
#'     )
#' )
#'
#' ## Summarise using summarise_abi_file()
#' processed_seq <- secondary_peaks(s4_sangerseq)
#'
#' @export
secondary_peaks <- function(
    s, ratio = 0.33,
    output.folder = NA,
    file.prefix = "seq",
    processors = NULL) {
    basecalls <- makeBaseCalls(s, ratio = ratio)

    primary <- primarySeq(basecalls, string = TRUE)
    secondary <- secondarySeq(basecalls, string = TRUE)


    comp <- compareStrings(primary, secondary)
    diffs <- str_locate_all(pattern = "\\?", comp)[[1]][, 1]
    primary.vector <- strsplit(primary, split = "")[[1]]
    secondary.vector <- strsplit(secondary, split = "")[[1]]

    primary.basecall <- primary.vector[diffs]
    secondary.basecall <- secondary.vector[diffs]

    r <- data.frame(
        "position" = diffs,
        "primary.basecall" = primary.basecall,
        "secondary.basecall" = secondary.basecall
    )

    if (!is.na(output.folder)) {
        if (dir.exists(output.folder)) {
            chromname <- paste(file.prefix, "_", "chromatogram.pdf", sep = "")
            chrom <- chromatogram(basecalls,
                width = 50, height = 2, showcalls = "both",
                filename = file.path(output.folder, chromname),
                trim5 = 100,
                trim3 = nchar(primary) - 150,
            )
        } else {
            warning(sprintf(
                "Couldn't find directory '%s', no files saved",
                output.folder
            ))
        }
    }

    return(list("secondary.peaks" = r, "read" = basecalls))
}

trim.mott <- function(abif.seq, cutoff = 0.0001) {
    if (!is(cutoff, "numeric") | cutoff < 0) {
        stop("cutoff must be a number of at least 0")
    }
    if (!is(abif.seq, "abif")) {
        stop("abif.seq must be an 'abif'
             object from the sangerseqR package")
    }
    abif.seq <- abif.seq@data
    start <- FALSE ## flag for starting position of trimmed sequence
    trim_start <- 0 ## init start index
    seqlen <- nchar(abif.seq$PBAS.2)
    qual <- abif.seq$PCON.2
    ## calculate base score
    score_list <- cutoff - (10**(qual / -10.0))
    ## calculate cummulative score
    ## if cumulative value < 0, set it to 0
    score <- score_list[1]
    if (score < 0) {
        score <- 0
    } else {
        trim_start <- 1
        start <- TRUE
    }
    cummul_score <- c(score)
    for (i in seq_along(score_list)[-1]) {
        score <- cummul_score[length(cummul_score)] + score_list[i]
        if (score <= 0) {
            cummul_score <- c(cummul_score, 0)
        } else {
            cummul_score <- c(cummul_score, score)
            if (start == FALSE) {
                trim_start <- i
                start <- TRUE
            }
        }
        ## trim_finish=index of highest cummulative score,
        ## marking the end of sequence segment with highest cummulative score
        trim_finish <- which.max(cummul_score)
    }

    ## fix an edge case, where all scores are worse than the cutoff
    ## in this case you wouldn't want to keep any bases at all
    if (sum(cummul_score) == 0) {
        trim_finish <- 0
    }
    return(list("start" = trim_start, "finish" = trim_finish))
}

fix.trims <- function(trims, seq.sanger, seq.abif, processors) {
    if (trims$start == 0 & trims$finish == 0) {
        return(trims)
    }

    ## First we trim the original sequence
    original.seq <- seq.abif@data$PBAS.2

    original.trimmed <- substring(original.seq, trims$start, trims$finish)

    ## Align the original and recalled sequences
    recalled <- primarySeq(seq.sanger, string = TRUE)
    seqs <- DNAStringSet(c(original.trimmed, recalled))
    pa <- AlignSeqs(seqs,
        iterations = 0, refinements = 0,
        verbose = FALSE, processors = processors
    )

    ## Get the sequence out, and find the first and last gaps.
    aligned.trimmed <- as.character(pa[[1]])
    not.gaps <- str_locate_all(aligned.trimmed, pattern = "[^-]")[[1]][, 1]

    start <- min(not.gaps)
    finish <- max(not.gaps)

    if (start < 1) {
        start <- 1
    }
    if (finish > nchar(recalled)) {
        finish <- nchar(recalled)
    }

    return(list("start" = start, "finish" = finish))
}
