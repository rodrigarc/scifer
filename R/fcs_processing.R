#' Extract index sorting information from flow cytometry data
#'
#' @param folder_path Folder containing all the flow data index filex (.fcs). Files should be named with their sample/plate ID.   eg. "E11_01.fcs"
#' @param compensation Logical argument, TRUE or FALSE, to indicate if the index files were compensated or not. If TRUE, it will apply its compensation prior assigning specificities
#' @param plate_wells Type of plate used for single-cell sorting. eg. "96" or "384"
#' @param probe1 Name of the first channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
#' @param probe2 Name of the second channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
#' @param posvalue_probe1 Threshold used for fluorescence intensities to be considered as positive for the first probe
#' @param posvalue_probe2 Threshold used for fluorescence intensities to be considered as positive for the second probe
#'
#' @return If saved as an object, it returns a table containing all the processed flow cytometry index files, with their fluorescence intensities for each channel and well position.
#' At the same time, it also plots a traditional flow density plot with the sorted cells and the selected thresholds for the two probes.
#'
#' @import dplyr ggplot2
#' @rawNamespace import(flowCore, except = filter)
#' @importFrom data.table rbindlist
#' @importFrom plyr mapvalues
#' @importFrom stats na.omit
#' @importFrom rlang .data
#' @importFrom scales trans_breaks trans_format math_format
#'
#' @examples
#' \dontrun{
#' index_sort_data <- fsc.processing(folder_path = "test/test_dataset/fcs_files/",
#' compensation = TRUE, plate_wells = 96,
#' probe1 = "Pre.F", probe2 = "Post.F",
#' posvalue_probe1 = 600, posvalue_probe2 = 400)
#'}
#' @export
fcs_processing <- function(folder_path = "test/test_dataset/fcs_files/",
                           compensation = TRUE, plate_wells = 96,
                           probe1 = "Pre.F", probe2 = "Post.F",
                           posvalue_probe1 = 600, posvalue_probe2 = 400) {
  fs <- read.flowSet(path = folder_path, truncate_max_range = FALSE)

  if (compensation == TRUE) {
    comp <- fsApply(fs, function(x) spillover(x)[[1]], simplify = FALSE)
    fs_comp <- compensate(fs, comp)
    print("Samples were compensated using the compensation saved on fsc index files.")
  } else if (compensation == FALSE) {
    fs_comp <- fs
    print("Samples were not compensated.")
  } else {
    print("Compensation argument should be TRUE or FALSE.")
  }

  fs_comp <- fsApply(fs_comp, simplify = FALSE, function(x) {
    # Change name of channels for marker's name
    params <- parameters(x)[["desc"]]
    colnames(x)[!is.na(params)] <- na.omit(params)

    # Extract single-cell sorted plate position
    if (plate_wells == 96) {
      df_fs_comp <- getIndexSort(x) %>% mutate(
        row = plyr::mapvalues(.data$XLoc, from = seq(0, 7), to = LETTERS[1:8]),
        column = plyr::mapvalues(.data$YLoc, from = seq(0, 11), to = sprintf("%02d", as.numeric(seq(1:12)))),
        well_ID = paste0(.data$row, .data$column)
      )
      print("96-well plates were used for sorting.")
    } else if (plate_wells == 384) {
      df_fs_comp <- getIndexSort(x) %>% mutate(
        row = plyr::mapvalues(.data$XLoc, from = seq(0, 15), to = LETTERS[1:16]),
        column = plyr::mapvalues(.data$YLoc, from = seq(0, 23), to = sprintf("%02d", as.numeric(seq(1:24)))),
        well_ID = paste0(.data$row, .data$column)
      )
      print("384-well plates were used for sorting.")
    } else {
      print("Only 96 or 384-well plates are supported")
    }

    return(df_fs_comp)
  })

  joined_fsc_table <- data.table::rbindlist(fs_comp, idcol = TRUE) %>%
    rename(sample_ID = .data$.id) %>%
    mutate(sample_ID = gsub(".fcs", "", .data$sample_ID))

  # Classify sorted single-cells according to specificity
  joined_fsc_table <- joined_fsc_table %>%
    mutate(specificity = case_when(
      get(probe1) > posvalue_probe1 & get(probe2) > posvalue_probe2 ~ "DP",
      get(probe1) > posvalue_probe1 & get(probe2) < posvalue_probe2 ~ probe1,
      get(probe1) < posvalue_probe1 ~ probe2
    ))
  # Plot classification according to selected thresholds
  .x <- NULL
  plot_thresholds <- joined_fsc_table %>%
    ggplot(aes(x = get(probe1), y = get(probe2))) +
    geom_point(size = .7, color = "grey40") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    geom_vline(xintercept = posvalue_probe1) +
    geom_hline(yintercept = posvalue_probe2) +
    labs(x = probe1, y = probe2) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(colour = "black", size = 1),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12)
    ) +
    geom_density2d(color = "black") +
    annotation_logticks()
  print(plot_thresholds)
  return(joined_fsc_table)
}
