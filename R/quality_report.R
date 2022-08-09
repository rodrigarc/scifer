#' Generate general and individualized reports
#'
#' This function uses the other functions already described to create a HTML report based on sequency quality. Besides the HTML reports, it also creates fasta files with all the sequences and individualized sequences, in addition to a csv file with the quality scores and sequences considered as good quality.
#'
#' @param data_folder Full file directory for searching all ab1 files in a recursive search method (all files including files in subfolders)
#' @param outputfile Output file name for the report generation
#' @param output_dir Output directory for all the different output files that are generated during the report
#' @param processors Number of processors to use, or NULL (the default) for all available processors
#' @param plot_chromatogram Logical argument, TRUE or FALSE, to indicate if chromatograms should be plotted or not. Default is FALSE.
#' @param folder_path_fcs  Full file directory for searching all flow cytometry index files (.fcs) in a recursive search method (all files including files in subfolders)
#' @return HTML reports, fasta files, csv file
#'
#' @examples
#' \dontrun{
#' quality.report(
#'   data_folder = "test/test_dataset/sanger_sequences",
#'   outputfile = "QC-report.html",
#'   output_dir = "test/",
#'   processors = 1
#' )
#' }
#' @importFrom rmarkdown render
#' @importFrom sangerseqR readsangerseq
#' @importFrom gridExtra grid.arrange
#' @importFrom stringr str_detect
#' @importFrom kableExtra kable kable_styling
#' @importFrom tibble rownames_to_column
#' @importFrom gridExtra grid.arrange
#' @importFrom sangerseqR primarySeq
#' @import ggplot2
#'
#' @export
quality_report <- function(data_folder, outputfile = "QC_report.html", output_dir = "test/", processors = 1, folder_path_fcs = "path/to/fcs_datasets",plot_chromatogram = FALSE) {
  input <- system.file("rmd", "HC_report.Rmd", package = "scifer")

  render(input,
    output_dir = output_dir,
    params = list(
      data_folder = data_folder,
      output_dir = output_dir,
      processors = processors,
      plot_chromatogram = plot_chromatogram,
      folder_path_fcs = folder_path_fcs
    ),
    output_file = outputfile
  )
}
