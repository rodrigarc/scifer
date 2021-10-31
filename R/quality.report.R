#' Generate general and individualized reports
#'
#' This function uses the other functions already described to create a HTML report based on sequency quality. Besides the HTML reports, it also creates fasta files with all the sequences and individualized sequences, in addition to a csv file with the quality scores and sequences considered as good quality.
#'
#' @param data_folder Full file directory for searching all ab1 files in a recursive search method (all files including files in subfolders)
#' @param outputfile Output file name for the report generation
#' @param output_dir Output directory for all the different output files that are generated during the report
#' @return HTML reports, fasta files, csv file
#'
## Examples not added due to error regarding location problems
## @examples
## quality.report(data_folder = "test/teste_dataset/",
##              outputfile = paste0(Sys.Date(),"_HC-QC-report.html"),
##              output_dir = "test/")
#'
#' @export quality.report

quality.report <- function(data_folder, outputfile, output_dir, processors = 1){
  input <- system.file("rmd", "HC_report.Rmd", package = "RepertoiR")

  rmarkdown::render(input,
                    output_dir = output_dir,
                    params = list(
                      data_folder = data_folder,
                      output_dir = output_dir,
                      processors = processors
                    ),
                    output_file = outputfile)
}

