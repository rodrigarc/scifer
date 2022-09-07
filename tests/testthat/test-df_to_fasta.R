test_that("output BStringSet", {
  suppressMessages(
    string <- df_to_fasta(sequence_name=c("myseq1"),
                             sequence_strings=c("AATGTCTG"),save_fasta = FALSE)
  )
  expect_s4_class(string, "BStringSet")
})

test_that("different number of names and sequences", {
  suppressMessages(
    expect_warning(df_to_fasta(sequence_name=c("myseq1"),
                               sequence_strings=c("AATGTCTG","ATGTCA"),save_fasta = FALSE),
                   regexp = "Sequence column does not have the same length as sequences name")
  )
})

test_that("not save fasta files", {
  expect_message(df_to_fasta(sequence_name=c("myseq1"),
                             sequence_strings=c("AATGTCTG"),save_fasta = FALSE),
                 regexp = "Fasta file not saved.")
})

test_that("save fasta files", {

  t <- tempdir()
  fs::dir_ls(t) %>% fs::file_delete()
  df_to_fasta(sequence_name=c("myseq1"),
              sequence_strings=c("AATGTCTG"),
              output_dir = t,
              file_name = "fasta_test.fasta",
              save_fasta = TRUE)
  expect_true(fs::dir_ls(t) %>% length() > 0)
  fs::dir_ls(t) %>% fs::file_delete()
})
