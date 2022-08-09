#' Fasta file creation from dataframe columns and/or vectors.
#'
#' @param sequence_name Vector containing the names for each sequence, usually a column from a data.frame. eg. df$sequence_name
#' @param sequence_strings Vector containing the DNA or RNA or AA sequences, usually a column from a data.frame. eg. df$sequences
#' @param file_name Output file name to be saved as a fasta file
#' @param output_dir Output directory for the fasta file. Default is the working directory
#'
#' @return Saves a fasta file in the desired location, and also returns the stringset as BStringSet if saved as an object.
#' @importFrom Biostrings BStringSet writeXStringSet
#' @export
#' @examples
#' \dontrun{
#' df_to_fasta(sequence_name = c("myseq1", "myseq2"),
#' sequence_strings = c("GATCGAT","ATCGTAG"),
#' file_name = "my_fasta.fasta",
#' output_dir = "path/to/output")
#'
#' df_to_fasta(sequence_name = df$sequence_name,
#' sequence_strings = df$sequence,
#' file_name = "my_fasta.fasta",
#' output_dir = "path/to/output")
#' }
df_to_fasta <- function(sequence_name, sequence_strings, file_name = "sequences.fasta", output_dir = NULL) {
  if (length(sequence_strings) != length(sequence_name)) {
    print("Sequences columng does not have the same length as sequences name")
  } else {
    str <- BStringSet(sequence_strings)
    names(str) <- sequence_name
  }
  if (is.null(output_dir)) {
    writeXStringSet(str, filepath = file_name, append = FALSE, format = "fasta")
  } else {
    writeXStringSet(str, filepath = paste0(output_dir, "/", file_name), append = FALSE, format = "fasta")
  }
  return(str)
}
