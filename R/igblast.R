#' Run IgDiscover for IgBlast using basilisk, which enables the python environment for Igblast
#'
#' @param database Vector containing the database for VDJ sequences
#' @param fasta Vector containing the sequences, usually a column from a data.frame. eg. df$sequences
#' @param threads Variable containing the number of cores when computing in parallel, default threads = 1
#' 
#' @return Creates a data frame with the Igblast analysis where each row is the tested sequence with columns containing the results for each sequence
#' @import reticulate here basilisk 
#' @importFrom utils read.table
#' @examples
#' ## Example with test sequences
#' igblast(
#'     database = system.file("/inst/extdata/test_fasta/KIMDB_rm", package = "scifer"),
#'     fasta = system.file("/inst/extdata/test_fasta/test_igblast.txt", package = "scifer"),
#'     threads = 1
#' )
#' @export
igblast <- function(database="path/to/folder", fasta="path/to/file", threads = 1) {
    db_dir <- here(database)
    fasta_dir <- here(fasta) 
    if (!dir.exists(db_dir) | is.null(db_dir) | is.na(db_dir) | db_dir == "" | length(db_dir) > 1) {
      stop(paste("The database directory does not exist."))
    }
    if (!file.exists(fasta_dir) | is.null(fasta_dir) | is.na(fasta_dir) | fasta_dir == "" | length(fasta_dir) > 1) {
      stop(paste("The fasta file directory", fasta_dir, "does not exist."))
    } 
    proc <- basiliskStart(env)
    on.exit(basiliskStop(proc))
    py_script <- system.file("script/igblastwrap.py", package = "scifer")
      run_igblastwrap <- basiliskRun(proc, fun=function(arg1, arg2, arg3) {
        tryCatch({
          stderr_file <- tempfile()
          sink(stderr_file, type="message")
          try(source_python(py_script), silent=TRUE)
            } , error = function(e) NULL 
          )
        df <- system(paste("python", py_script, "--threads", threads, database, fasta,  sep = " "),
                    intern = TRUE)
        if (length(df) == 0 || is.null(df) || all(is.na(df)) || all(df == ""))  {
            warning("Data frame is empty. Sequences not aligned.")
            results_airr <- NULL  
          } else {
            results_airr <- utils::read.table(text=df, header=TRUE, sep = "\t")
          }  
        return(final_output = results_airr)
      }, arg1=database, arg2=fasta, arg3=threads)
    run_igblastwrap
}