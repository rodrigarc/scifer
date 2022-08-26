test_that("output BStringSet", {
  string <- df_to_fasta(sequence_name=c("myseq1"),
                           sequence_strings=c("AATGTCTG"),save_fasta = FALSE)
  expect_condition(string, class = "BStringSet")

})
