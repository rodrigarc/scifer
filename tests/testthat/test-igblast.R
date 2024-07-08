test_that("the database directory does not exists NA", {
    expect_error(
        igblast(database = NA, fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
        "The database directory does not exist.")
})

test_that("the database directory does not exists INVALID", {
    expect_error(
        igblast(database = "/invalidpath", fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
        "The database directory does not exist.")
})

test_that("the database directory does not exists SPACE", {
    expect_message(igblast(database = "", fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
        "Data frame is empty. Sequences not aligned.") 
})

test_that("the database directory does not exists NULL", {
    expect_error(
        igblast(database = NULL, fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
        "argument is of length zero")
})

test_that("the file directory does not exists NA", {
    expect_error(
        igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), fasta = NA, threads = 1),
        "The fasta file directory does not exist.")
})

test_that("the file directory does not exists INVALID", {
    expect_error(
        igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), fasta = "/invalidpath", threads = 1),
        "The fasta file directory does not exist.")
})

test_that("the file directory does not exists SPACE", {
    expect_message(
        igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), fasta = "", threads = 1),
        "Data frame is empty. Sequences not aligned.")
})

test_that("the file directory does not exists NULL", {
    expect_error(
        igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), fasta = NULL, threads = 1),
        "argument is of length zero")
})

test_that("returns a data.frame object", {
    result <- igblast(
                    database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), 
                    system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), 
                    threads = 1
                    )
    expect_s3_class(result, "data.frame")
})
