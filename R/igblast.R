#' Run IgBLAST python wrapper
#'
#' A wrapper funtion to run the IgBLAST python script to annotate VDJ sequences.
#' It is python based and relies on conda environments that are built when the
#' funtion is called.
#' 
#' @param database Vector containing the database for VDJ sequences
#' @param fasta Vector containing the sequences, usually a column from a data.frame. eg. df$sequences
#' @param threads Variable containing the number of cores when computing in parallel, default threads = 1
#'
#' @return Creates a data frame with the IgBLAST annotation where each row is the queried sequence with columns containing the IgBLAST results
#' @import reticulate here basilisk basilisk.utils
#' @importFrom utils read.table
#' @examples
#' ## Example with test sequences
#' \dontrun{
#' igblast(
#'     database = system.file("/extdata/test_fasta/KIMDB_rm", package = "scifer"),
#'     fasta = system.file("/extdata/test_fasta/test_igblast.txt", package = "scifer"),
#'     threads = 1
#' )
#' }
#' @export
igblast <- function(database, fasta, threads = 1) {
  # Input validation BEFORE path resolution
  if (is.null(database) || length(database) == 0 || is.na(database) || database == "")
    stop("The database directory does not exist.")
  if (is.null(fasta) || length(fasta) == 0 || is.na(fasta) || fasta == "")
    stop("The fasta file directory does not exist.")
  if (!is.numeric(threads) || length(threads) != 1)
    stop("The threads argument should be a numeric value.")
  
  db <- here::here(database)
  fa <- here::here(fasta)
  
  if (!dir.exists(db)) stop("The database directory does not exist.")
  if (!file.exists(fa)) stop("The fasta file directory does not exist.")
  # Environment specification
  env_spec <- list(
    packages = c(
      "python==3.9.19", "igblast==1.22.0", "cffi==1.16.0", "python-isal==1.6.1",
      "pip==24.0", "pycparser==2.22", "setuptools==70.1.1", "wheel==0.43.0",
      "xopen==2.0.1", "python-zlib-ng==0.4.3", "zstandard==0.22.0", "dnaio==1.2.1"
    ),
    channels = c("bioconda", "conda-forge")
  )
  
  # Check OS support
  if (!(isLinux() || (isMacOSX() && !isMacOSXArm()))) {
    message("igblast environment setup not supported on this OS/architecture. Skipping.")
    return(invisible(NULL))
  }
  
  # Create environment if needed
  env_name <- "igblast_wrap_basilisk"
  env_path <- basilisk.utils::createEnvironment(
    pkg = "scifer",
    name = env_name,
    version = "1.0.0",
    packages = env_spec$packages,
    channels = env_spec$channels
  )
  
  # Path to python script
  py_script <- system.file("script/igblastwrap.py", package = "scifer")
  
  # Run python script inside the created conda env using 'conda run'
  cmd <- c("run", "--prefix", env_path, "python", py_script,
           "--threads", as.character(threads),
           "--database", db,
           "--fasta", fa)
  
  res <- tryCatch(
    system2("conda", cmd, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      message("Failed to run igblast: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(res) || length(res) == 0 || all(res == "")) {
    message("No output returned from igblast.")
    return(NULL)
  }
  
  # Parse output to data.frame
  df <- tryCatch({
    header_line <- grep("^sequence_id", res)
    txt <- paste(res[header_line:length(res)][res[header_line:length(res)] != ""], collapse = "\n")
    read.table(text = txt, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    message("Failed to parse igblast output: ", e$message)
    NULL
  })
  
  return(df)
}
