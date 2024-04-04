#' Extract index sorting information from flow cytometry data
#'
#' @param folder_path Folder containing all the flow data index filex (.fcs). Files should be named with their sample/plate ID.   eg. "E11_01.fcs"
#' @param compensation Logical argument, TRUE or FALSE, to indicate if the index files were compensated or not. If TRUE, it will apply its compensation prior assigning specificity
#' @param plate_wells Type of plate used for single-cell sorting. eg. "96" or "384"
#' @param probe1 Name of the first channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
#' @param probe2 Name of the second channel used for the probe or the custom name assigned to the channel in the index file. eg. "FSC.A", "FSC.H", "SSC.A","DsRed.A", "PE.Cy5_5.A", "PE.Cy7.A","BV650.A", "BV711.A","Alexa.Fluor.700.A" "APC.Cy7.A","PerCP.Cy5.5.A","Time"
#' @param posvalue_probe1 Threshold used for fluorescence intensities to be considered as positive for the first probe
#' @param posvalue_probe2 Threshold used for fluorescence intensities to be considered as positive for the second probe
#'
#' @return If saved as an object, it returns a table containing all the processed flow cytometry index files, with their fluorescence intensities for each channel and well position.
#'
#' @import dplyr
#' @rawNamespace import(flowCore, except=c(filter, cytolib_LdFlags))
#' @importFrom data.table rbindlist
#' @importFrom plyr mapvalues
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @examples
#' index_sort_data <- fcs_processing(
#'     folder_path = system.file("/extdata/fcs_index_sorting",
#'                                 package = "scifer"),
#'     compensation = TRUE, plate_wells = 96,
#'     probe1 = "Pre.F", probe2 = "Post.F",
#'     posvalue_probe1 = 600, posvalue_probe2 = 400
#' )
#'
#' @export
fcs_processing <- function(folder_path = "test/test_dataset/fcs_files/",
    compensation = TRUE, plate_wells = 96,
    probe1 = "Pre.F", probe2 = "Post.F",
    posvalue_probe1 = 600, posvalue_probe2 = 400) {
    fs <- read.flowSet(path = folder_path, truncate_max_range = FALSE)
    if (compensation == TRUE) {
      if(all(unlist(lapply(spillover(fs[[1]]), is.null)))) {
        stop("Compensation matrices are all empty")
      } else {
        not_null_comp <- which(!unlist(lapply(spillover(fs[[1]]), is.null)))
        if (length(not_null_comp) > 1) {
          not_null_comp <- not_null_comp[1]
        }
        comp <- fsApply(fs, function(x) spillover(x)[[not_null_comp]],
                        simplify = FALSE)
        fs_comp <- compensate(fs, comp)
        message("Samples were compensated using the data on fsc index files.")
      }
    } else if (compensation == FALSE) {
        fs_comp <- fs
        message("Samples were not compensated.")
    } else {
        stop("Compensation argument should be TRUE or FALSE.")
    }
    fs_comp <- fsApply(fs_comp, simplify = FALSE, function(x) {
        ## Change name of channels for marker's name
        params <- parameters(x)[["desc"]]
        colnames(x)[!is.na(params)] <- na.omit(params)
        ## Extract single-cell sorted plate position
        if (plate_wells == 96) {
            df_fs_comp <- getIndexSort(x) %>% mutate(
                row = plyr::mapvalues(.data$XLoc,
                    from = seq(0, 7), to = LETTERS[seq_len(8)],
                    warn_missing = FALSE
                ),
                column = plyr::mapvalues(.data$YLoc,
                    from = seq(0, 11),
                    to = sprintf("%02d", as.numeric(seq_len(12))),
                    warn_missing = FALSE
                ),
                well_ID = paste0(.data$row, .data$column)
            )
            message("96-well plates were used for sorting.")
        } else if (plate_wells == 384) {
            df_fs_comp <- getIndexSort(x) %>% mutate(
                row = plyr::mapvalues(.data$XLoc,
                    from = seq(0, 15), to = LETTERS[seq_len(16)],
                    warn_missing = FALSE
                ),
                column = plyr::mapvalues(.data$YLoc,
                    from = seq(0, 23),
                    to = sprintf("%02d", as.numeric(seq_len(24))),
                    warn_missing = FALSE
                ),
                well_ID = paste0(.data$row, .data$column)
            )
            message("384-well plates were used for sorting.")
        } else {
            stop("Only 96 or 384-well plates are supported")
        }

        return(df_fs_comp)
    })
    joined_fsc_table <- data.table::rbindlist(fs_comp, idcol = TRUE) %>%
        rename(sample_ID = .data$.id) %>%
        mutate(sample_ID = gsub(".fcs", "", .data$sample_ID))
    ## Classify sorted single-cells according to specificity
    joined_fsc_table <- joined_fsc_table %>%
        mutate(specificity = case_when(
            get(probe1) < posvalue_probe1 &
                get(probe2) < posvalue_probe2 ~ "DN",
            get(probe1) > posvalue_probe1 &
                get(probe2) > posvalue_probe2 ~ "DP",
            get(probe1) > posvalue_probe1 &
                get(probe2) < posvalue_probe2 ~ probe1,
            get(probe1) < posvalue_probe1 &
                get(probe2) > posvalue_probe2 ~ probe2
        ))

    probes_values <- data.frame(
        probes = c(probe1, probe2),
        posvalues = c(posvalue_probe1, posvalue_probe2)
    )

    return(list("processed_fcs" = joined_fsc_table,
                "selected_probes" = probes_values))
}
