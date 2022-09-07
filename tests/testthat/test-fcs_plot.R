test_that("Output ggplot", {
  suppressMessages(
  index_sort <- fcs_processing(
    folder_path=system.file("/extdata/fcs_index_sorting",package = "scifer"),
    compensation=TRUE, plate_wells = 96)
  )
  suppressWarnings(
    g1 <- fcs_plot(index_sort)
    )
  expect_s3_class(g1, "ggplot")
})

test_that("Input list from fcs_processing()", {
  expect_error(fcs_plot(), "Input data.frame is NULL/empty")
  expect_error(fcs_plot(2), "Input is not a list object")
  expect_error(fcs_plot(list("1")), "Input is not a list with length 2")
  expect_error(fcs_plot(list("1","2","3")), "Input is not a list with length 2")
  expect_error(fcs_plot(list(list("1"),"2")), "At least one of the objects inside the list is not a data.frame")
  expect_error(fcs_plot(list(as.matrix("1"),"2")), "At least one of the objects inside the list is not a data.frame")
})
