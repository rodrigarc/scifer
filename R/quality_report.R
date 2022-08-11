#' Generate general and individualized reports
#'
#' This function uses the other functions already described to create a HTML report based on sequencing quality. Besides the HTML reports, it also creates fasta files with all the sequences and individualized sequences, in addition to a csv file with the quality scores and sequences considered as good quality.
#'
#' @param folder_sequences Full file directory for searching all ab1 files in a recursive search method. It includes all files in subfolders
#' @param outputfile Output file name for the report generation
#' @param output_dir Output directory for all the different output files that are generated during the report
#' @param processors Number of processors to use, you can set to NULL to detect automatically all available processors
#' @param plot_chromatogram Logical argument, TRUE or FALSE, to indicate if chromatograms should be plotted or not. Default is FALSE
#' @param folder_path_fcs  Full file directory for searching all flow cytometry index files, files with .fcs extensions, in a recursive search method
#' @param raw_length Minimum sequence length for filtering. Default is 400 for B cell receptors
#' @param trim_start Starting position where the sequence should start to have a good base call accuracy. Default is 50 for B cell receptors
#' @param trim_finish Last position where the sequence should have a good base call accuracy. Default is 409 for B cell receptors
#' @param trimmed_mean_quality Minimum Phred quality score expected for an average sequence. Default is 30, which means average of 99.9\% base call accuracy
#' @param compensation Logical argument, TRUE or FALSE, to indicate if the index files were compensated or not. If TRUE, it will apply its compensation prior assigning specificities
#' @param plate_wells Type of plate used for single-cell sorting. eg. "96" or "384"
#' @param probe1 Name of the first channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
#' @param probe2 Name of the second channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
#' @param posvalue_probe1 Threshold used for fluorescence intensities to be considered as positive for the first probe
#' @param posvalue_probe2 Threshold used for fluorescence intensities to be considered as positive for the second probe
#' @param cdr3_start Expected CDR3 starting position, that depends on your primer set. Default is position 100
#' @param cdr3_end Expected CDR3 end position, that depends on your primer set. Default is position 150
#'
#'
#' @return Saves HTML reports, fasta files, csv files
#'
#' @examples
#' \dontrun{
#' quality_report(
#'     data_folder="~/test/test_dataset/sanger_sequences",
#'     outputfile="QC-report.html",
#'     output_dir="~/test/",
#'     folder_path_fcs="~/path/to/fcs_datasets",
#'     processors=4, compensation=TRUE, plate_wells="96",
#'     probe1="Pre.F", probe2="Post.F",
#'     posvalue_probe1=600, posvalue_probe2=400,
#'     cdr3_start=100,
#'     cdr3_end=150
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
#' @import ggplot2 dplyr knitr
#'
#' @examples
#' \dontrun{
#' quality_report(
#' folder_sequences=system.file("extdata/sorted_sangerseq/", package="scifer"),
#' output_dir="", folder_path_fcs=system.file("extdata/fcs_index_sorting/", package="scifer"),
#'  probe1="Pre.F", probe2="Post.F", posvalue_probe1=600, posvalue_probe2=400
#' )}
#'
#' @export
quality_report <- function(folder_sequences = "path/to/sanger_sequences", outputfile="QC_report.html", output_dir="test/", processors=NULL,
    folder_path_fcs="path/to/fcs_datasets", plot_chromatogram=FALSE,
    raw_length=400, trim_start=50, trim_finish=409,
    trimmed_mean_quality=30,
    compensation=TRUE, plate_wells="96",
    probe1="Pre.F", probe2="Post.F",
    posvalue_probe1=600, posvalue_probe2=400,
    cdr3_start=100,
    cdr3_end=150) {

    input <- system.file("rmd", "HC_report.Rmd", package="scifer")

    render(input,
        output_dir=output_dir,
        params=list(
            folder_sequences=folder_sequences,
            output_dir=output_dir,
            processors=processors,
            plot_chromatogram=plot_chromatogram,
            folder_path_fcs=folder_path_fcs,
            raw_length=raw_length,
            trim_start=trim_start,
            trim_finish=trim_finish,
            trimmed_mean_quality=trimmed_mean_quality,
            compensation=compensation,
            plate_wells=plate_wells,
            probe1=probe1,
            probe2=probe2,
            posvalue_probe1=posvalue_probe1,
            posvalue_probe2=posvalue_probe2,
            cdr3_start=cdr3_start,
            cdr3_end=cdr3_end
        ),
        output_file=outputfile
    )
}
