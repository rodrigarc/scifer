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
#'                                 
#'     fasta = system.file("/inst/extdata/test_fasta/test_igblast.txt", package = "scifer"),
#'     
#'     threads = 1
#' )
#' @export

igblast <- function(database = "path/to/folder", fasta = "path/to/file", threads = 1) {
  db_dir <- here(database)
  fasta_dir <- here(fasta) 

  #db_dir <- here("inst", "extdata", "test_fasta", "KIMDB_rm")
  #fasta_dir <- here("inst", "extdata", "test_fasta")
  
  if (dir.exists(db_dir)) {
    #print(paste("The database directory", db_dir, "exists."))
  } else {
    stop(paste("The database directory", db_dir, "does not exist."))
  }
  
  if (file.exists(fasta_dir)) {
   # print(paste("The fasta file directory", fasta_dir, "exists."))
  } else {
    stop(paste("The fasta file directory", fasta_dir, "does not exist."))
  }
  
  proc <- basiliskStart(env)
  on.exit(basiliskStop(proc))
  
  py_script <- system.file("script/igblastwrap.py", package = "scifer")
  
  run_igblastwrap <- basiliskRun(proc, fun=function(arg1, arg2, arg3) {
    
   try(source_python(py_script), silent = TRUE)
   #try(reticulate::source_python("/script/igblastwrap.py"), silent = TRUE)
   #source_python(py_script)
   #source_python("/script/igblastwrap.py")
      df <- system(paste("python", py_script, "--threads", threads, database, fasta,  sep = " "),
                  intern = TRUE)
#if + warning if empty 
      if (length(df) == 0) {
        warning("The data frame is empty. Sequences were not aligned.")
        results_airr <- NULL  # or handle appropriately, e.g., create an empty data frame
      } else {
        results_airr <- utils::read.table(text = df, header = TRUE, sep = "\t")
      }  
      
    #results_airr <- utils::read.table(text = df, header = TRUE, sep = "\t")

    return(final_output = results_airr)
      #}
  }, arg1=database, arg2=fasta, arg3=threads)
  
  run_igblastwrap
}



