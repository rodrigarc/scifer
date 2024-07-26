#' Run IgDiscover for IgBlast using basilisk, which enables the python environment for Igblast
#'
#' @param database Vector containing the database for VDJ sequences
#' @param fasta Vector containing the sequences, usually a column from a data.frame. eg. df$sequences
#' @param threads Variable containing the number of cores when computing in parallel, default threads = 1
#'
#' @return Creates a data frame with the Igblast analysis where each row is the tested sequence with columns containing the results for each sequence
#' @import reticulate here basilisk basilisk.utils
#' @importFrom utils read.table
#' @examples
#' ## Example with test sequences
#' 
#' igblast(
#'     database = system.file("/extdata/test_fasta/KIMDB_rm", package = "scifer"),
#'     fasta = system.file("/extdata/test_fasta/test_igblast.txt", package = "scifer"),
#'     threads = 1
#' )
#' @export
igblast <- function(database = "path/to/folder", fasta = "path/to/file", threads = 1) {
    db_dir <- here(database)
    fasta_dir <- here(fasta)
    if (!dir.exists(db_dir) || is.null(db_dir) || is.na(db_dir) || length(db_dir) > 1 || length(db_dir) == 0  || database == "" ) {
        stop("The database directory does not exist.")
    } else if (!file.exists(fasta_dir) || is.null(fasta_dir) || is.na(fasta_dir) || length(fasta_dir) > 1 || length(fasta_dir) == 0 || fasta == "" ) {
        stop("The fasta file directory does not exist.")
    } else if (!is.numeric(threads) || is.null(threads) || is.na(threads) || length(threads) == 0) {
        stop("The threads argument should be a numeric value.")
    }
    
    if (basilisk.utils::isMacOSXArm() == TRUE) {
        conda_create("temp", packages = c("python=3.9"))
        system2("conda", args = c("congfig", "--env", "--set", "subdir", "osx-64"))
        use_condaenv("temp")
    }
    # conda_create("temp", packages = c("python=3.9"))
    # system2("conda", args = c("congfig", "--env", "--set", "subdir", "osx-64"))
    # use_condaenv("temp")
    
    proc <- basiliskStart(env)
    on.exit(basiliskStop(proc))
    py_script <- system.file("script/igblastwrap.py", package = "scifer")
    run_igblastwrap <- basiliskRun(proc, fun = function(arg1, arg2, arg3) {
        tryCatch(
            {
                stderr_file <- tempfile()
                sink(stderr_file, type = "message")
                try(source_python(py_script), silent = TRUE)
            },
            error = function(e) NULL
        )
        if(isWindows()){
          if(system("makeblastdb") %in% c(127, "Exit Code 127")){
            stop("IgBLAST is not available on this system. Please install IgBLAST.\nFor more information, refer to the FAQ section in scifer's GitHub README.")
          }
        }
        df <- system2("python",
                      args = c(py_script, "--threads", threads, database, fasta),
                      stdout = TRUE)
        if (length(df) == 0 || is.null(df) || all(is.na(df)) || all(df == "")) {
            message("Data frame is empty. Sequences not aligned.")
            results_airr <- NULL
        } else {
            if(basilisk.utils::isWindows()){
              results_airr <- utils::read.table(text = df[-1], header = TRUE, sep = "\t")
            } else {
              results_airr <- utils::read.table(text = df, header = TRUE, sep = "\t")
            }
        }
        return(final_output = results_airr)
    }, arg1 = database, arg2 = fasta, arg3 = threads)
    run_igblastwrap
}
