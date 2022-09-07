test_that("Output list with df", {
  abi_seq <- sangerseqR::read.abif(
    system.file("/extdata/sorted_sangerseq/E18_C1/A1_3_IgG_Inner.ab1", package="scifer"))
  output <- summarise_abi_file(abi_seq)

  expect_type(output, "list")
  expect_type(output[["summary"]], "double")
  expect_s3_class(output[["quality_score"]], "data.frame")
})

test_that("Input of sangerseq object", {
  expect_error(summarise_abi_file("not_a_file"))
  expect_error(summarise_abi_file(1))
})
