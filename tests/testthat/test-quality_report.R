test_that("Output files", {
  t <- tempdir()
  fs::dir_ls(t) %>% fs::file_delete()
  suppressWarnings(
    quality_report(
      folder_sequences=system.file("extdata/sorted_sangerseq/", package="scifer"),
      output_dir=t, folder_path_fcs=system.file("extdata/fcs_index_sorting/", package="scifer"),
      probe1="Pre.F", probe2="Post.F", posvalue_probe1=600, posvalue_probe2=400, processors = 1)
)
  # test if files were generated
  expect_true(fs::dir_ls(t) %>% length() > 0)
  fs::dir_ls(t) %>% fs::file_delete()
})
