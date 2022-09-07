#' Plot flow data from index sorted cells
#'
#' @param processed_fcs_list List generated using `fcs_processing()` containing two data.frames
#'
#' @return Returns a ggplot object with a traditional flow density plot with the sorted cells and the selected thresholds for the two probes used in fcs_processing().
#'
#' @import ggplot2
#' @importFrom scales trans_breaks trans_format math_format
#'
#' @examples
#' index_sort_data <- fcs_processing(
#'     folder_path=system.file("/extdata/fcs_index_sorting",package = "scifer"),
#'     compensation=TRUE, plate_wells=96,
#'     probe1="Pre.F", probe2="Post.F",
#'     posvalue_probe1=600, posvalue_probe2=400)
#'
#' fcs_plot_obj <- fcs_plot(index_sort_data)
#'
#' @export
fcs_plot <- function(processed_fcs_list=NULL) {
  if (is.null(processed_fcs_list)) {
    stop("Input data.frame is NULL/empty")
  } else if (!is(processed_fcs_list, "list")) {
    stop("Input is not a list object")
  } else if (length(processed_fcs_list) != 2) {
    stop("Input is not a list with length 2")
  } else if (any(!is(processed_fcs_list[[1]], "data.frame")|
                 !is(processed_fcs_list[[2]], "data.frame"))) {
    stop("At least one of the objects inside the list is not a data.frame")
  } else {
    ## Extract values for plotting
    probe1 <- processed_fcs_list[["selected_probes"]][1,1]
    probe2 <- processed_fcs_list[["selected_probes"]][2,1]
    posvalue_probe1 <- processed_fcs_list[["selected_probes"]][1,2]
    posvalue_probe2 <- processed_fcs_list[["selected_probes"]][2,2]
    ## Create object to remove NOTE error
    .x <- NULL
    flow_plot_obj <- ggplot(processed_fcs_list[["processed_fcs"]], aes(x=get(probe1), y=get(probe2))) +
      geom_point(size=.7, color="grey40") +
      scale_x_log10(
        breaks=scales::trans_breaks("log10", function(x) 10^x),
        labels=scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(
        breaks=scales::trans_breaks("log10", function(x) 10^x),
        labels=scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_vline(xintercept=posvalue_probe1) +
      geom_hline(yintercept=posvalue_probe2) +
      labs(x=probe1, y=probe2) +
      theme_bw() +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background=element_rect(colour="black", size=1),
            axis.text=element_text(size=10, color="black"),
            axis.title=element_text(size=12)) +
      geom_density2d(color="black") +
      annotation_logticks()
  }
 return(flow_plot_obj)
}
