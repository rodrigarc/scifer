library(RepertoiR)
library(flowCore)
library(scales)

fsc.processing <- function(file.path = "test/teste_dataset/fcs_file/",
                           compensation = TRUE, plate_wells = 96,
                           probe1 = "Pre.F", probe2 = "Post.F",
                           posvalue_probe1 = 600, posvalue_probe2 = 400) {
  fs <- read.flowSet(path = file.path, truncate_max_range = FALSE)

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
        row = plyr::mapvalues(XLoc, from = seq(0, 7), to = LETTERS[1:8]),
        column = plyr::mapvalues(YLoc, from = seq(0, 11), to = sprintf("%02d", as.numeric(seq(1:12)))),
        well_ID = paste0(row, column)
      )
      print("96-well plates were used for sorting.")
    } else if (plate_wells == 384) {
      df_fs_comp <- getIndexSort(x) %>% mutate(
        row = plyr::mapvalues(XLoc, from = seq(0, 15), to = LETTERS[1:16]),
        column = plyr::mapvalues(YLoc, from = seq(0, 23), to = sprintf("%02d", as.numeric(seq(1:24)))),
        well_ID = paste0(row, column)
      )
      print("384-well plates were used for sorting.")
    } else {
      print("Only 96 or 384-well plates are supported")
    }

    return(df_fs_comp)
  })

  joined_fsc_table <- data.table::rbindlist(fs_comp)

  # Classifying sorted single-cells according to specificity
  joined_fsc_table <- joined_fsc_table %>%
    mutate(specificity = case_when(
      get(probe1) > posvalue_probe1 & get(probe2) > posvalue_probe2 ~ "DP",
      get(probe1) > posvalue_probe1 & get(probe2) < posvalue_probe2 ~ probe1,
      get(probe1) < posvalue_probe1 ~ probe2
    ))

  plot_thresholds <- joined_fsc_table %>%
    ggplot(aes(x = get(probe1), y = get(probe2))) +
    geom_point(size = .7, color = "grey40") +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
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

x <- fsc.processing(
  file.path = "test/teste_dataset/fcs_file/",
  compensation = TRUE, plate_wells = 96,
  probe1 = "Pre.F", probe2 = "Post.F",
  posvalue_probe1 = 600, posvalue_probe2 = 400
)
