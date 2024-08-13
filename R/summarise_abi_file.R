#' Create a summary of a single ABI sequencing file
#'
#' Takes a single ABI sequencing file and returns a summary of the file. 
#' The summary includes basic quality control metric of the sequence.
#'
#' @param seq.abif an abif.seq s4 object from the sangerseqR package
#' @param trim.cutoff the cutoff at which you consider a base to be bad. This works on a logarithmic scale, such that if you want to consider a score of 10 as bad, you set cutoff to 0.1; for 20 set it at 0.01; for 30 set it at 0.001; for 40 set it at 0.0001; and so on. Contiguous runs of bases below this quality will be removed from the start and end of the sequence. Default is 0.0001.
#' @param secondary.peak.ratio the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not.
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' @param processors Number of processors to use, or NULL (the default) for all available processors
#'
#' @return A numeric vector including:
#'          \enumerate{
#'              \item {raw.length}: the length of the untrimmed sequence, note that this is the sequence after conversion to a sangerseq object, and then the recalling the bases with MakeBaseCalls from the sangerseqR package\cr
#'              \item {trimmed.length}: the length of the trimmed sequence, after trimming using trim.mott from this package and the parameter supplied to this function \cr
#'              \item {trim.start}: the start position of the good sequence, see trim.mott for more details\cr
#'              \item {trim.finish}: the finish position of the good sequence, see trim.mott for more details\cr
#'              \item {raw.secondary.peaks}: the number of secondary peaks in the raw sequence, called with the secondary.peaks function from this package and the parameters supplied to this function \cr
#'              \item {trimmed.secondary.peaks}: the number of secondary peaks in the trimmed sequence, called with the secondary.peaks function from this package and the parameters supplied to this function \cr
#'              \item {raw.mean.quality}: the mean quality score of the raw sequence \cr
#'              \item {trimmed.mean.quality}: the mean quality score of the trimmed sequence \cr
#'              \item {raw.min.quality}: the minimum quality score of the raw sequence \cr
#'              \item {trimmed.min.quality}: the minimum quality score of the trimmed sequence \cr
#'          }
#' @rawNamespace import(Biostrings, except=c(collapse, union, intersect, setdiff, setequal))
#' @import dplyr
#' @importFrom sangerseqR sangerseq primarySeq
#' @importFrom rlang .data
#'
#' @examples
#' ## Read abif using sangerseqR package
#' abi_seq <- sangerseqR::read.abif(
#'     system.file("/extdata/sorted_sangerseq/E18_C1/A1_3_IgG_Inner.ab1",
#'         package = "scifer"
#'     )
#' )
#'
#' ## Summarise using summarise_abi_file()
#' summarise_abi_file(abi_seq)
#'
#' @export
summarise_abi_file <- function(
    seq.abif,
    trim.cutoff = 0.0001,
    secondary.peak.ratio = 0.33,
    output.folder = NA,
    prefix = "seq",
    processors = NULL) {
    seq.sanger <- sangerseq(seq.abif)
    ## Get secondary peaks
    secondary.peaks.data <- scifer:::secondary_peaks(seq.sanger,
        secondary.peak.ratio,
        output.folder, prefix,
        processors = processors
    )
    secondary.peaks <- secondary.peaks.data[["secondary.peaks"]]
    seq.sanger <- secondary.peaks.data[["read"]]

    ## Trim sequence
    trims <- trim.mott(seq.abif, cutoff = trim.cutoff)
    qual <- seq.abif@data$PCON.2
    qual.trimmed <- qual[trims$start:trims$finish]

    if (length(qual.trimmed) == 0) {
        qual.trimmed <- c(NA)
    } ## Summarise

    ## Fix up the trim locations to correspond to sangerseq primaryseq object
    if (trims$start == 1 &&
        trims$finish == nchar(as.character(seq.abif@data$PBAS.2))) {
        ## Do nothing if sequence was not trimmed
        trim.start <- 1
        trim.finish <- length(primarySeq(seq.sanger))
    } else {
        trims.fixed <- fix.trims(trims, seq.sanger, seq.abif, processors)
        trim.start <- trims.fixed$start
        trim.finish <- trims.fixed$finish
    }

    ## Get trimmed and untrimmed version of raw data
    seq.trimmed <- seq.sanger@primarySeq[trim.start:trim.finish]
    secondary.peaks.trimmed <- secondary.peaks %>%
        filter(.data$position >= trim.start, .data$position <= trim.finish)

    read_summary <- c(
        "raw.length" = length(seq.sanger@primarySeq),
        "trimmed.length" = length(seq.trimmed),
        "trim.start" = trim.start,
        "trim.finish" = trim.finish,
        "raw.secondary.peaks" = nrow(secondary.peaks),
        "trimmed.secondary.peaks" = nrow(secondary.peaks.trimmed),
        "raw.mean.quality" = mean(qual),
        "trimmed.mean.quality" = mean(qual.trimmed),
        "raw.min.quality" = min(qual),
        "trimmed.min.quality" = min(qual.trimmed)
    )
    qual_position <- cbind.data.frame(
        "score" = as.numeric(seq.abif@data$PCON.2),
        "position" = c(seq_len(length(seq.abif@data$PCON.2)))
    )
    return(list("summary" = read_summary, "quality_score" = qual_position))
}
