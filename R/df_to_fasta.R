#' Fasta file creation from dataframe columns and/or vectors.
#'
#' Creates a fasta file from vectors of names and sequences.
#' 
#' @param sequence_name Vector containing the names for each sequence, usually a column from a data.frame. eg. df$sequence_name
#' @param sequence_strings Vector containing the DNA or RNA or AA sequences, usually a column from a data.frame. eg. df$sequences
#' @param file_name Output file name to be saved as a fasta file
#' @param output_dir Output directory for the fasta file. Default is the working directory
#' @param save_fasta Logical argument, TRUE or FALSE, to indicate if fasta files should be saved. Default is TRUE.
#'
#' @return Saves a fasta file in the desired location, and also returns the stringset as BStringSet if saved as an object.
#' @importFrom Biostrings BStringSet writeXStringSet
#' @examples
#' ## Example with vectors, default for save_fasta ir TRUE
#' df_to_fasta(
#'     sequence_name = c("myseq1", "myseq2"),
#'     sequence_strings = c("GATCGAT", "ATCGTAG"),
#'     file_name = "my_sequences.fasta",
#'     output_dir = "",
#'     save_fasta = FALSE
#' )
#' @export
df_to_fasta <- function(
    sequence_name,
    sequence_strings,
    file_name = "sequences.fasta",
    output_dir = NULL, save_fasta = TRUE) {
    if (length(sequence_strings) != length(sequence_name)) {
        warning("Sequence column has different length of sequences name")
    } else {
        str <- BStringSet(sequence_strings)
        names(str) <- sequence_name
    }

    if (isTRUE(save_fasta)) {
        if (is.null(output_dir)) {
            writeXStringSet(str,
                filepath = file_name,
                append = FALSE, format = "fasta"
            )
        } else {
            writeXStringSet(str,
                filepath = paste(output_dir, file_name, sep = "/"),
                append = FALSE, format = "fasta"
            )
        }
    } else {
        message("Fasta file not saved.")
    }
    return(str)
}
