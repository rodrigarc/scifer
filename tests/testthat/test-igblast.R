# test_that("the database directory does not exists NULL", {
#   expect_error(igblast(database = NULL, fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
#                "The database directory does not exist.")
# })

test_that("the database directory does not exists NA", {
    expect_error(
        igblast(database = NA, fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
        "The database directory does not exist."
    )
})

test_that("the database directory does not exists INVALID", {
    expect_error(
        igblast(database = "/invalidpath", fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
        "The database directory does not exist."
    )
})

# test_that("the database directory does not exists space", {
#   expect_error(igblast(database = "", fasta = system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1),
#                "The database directory does not exist.")
# })

test_that("the file directory does not exists", {
    expect_error(
        igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), fasta = NA, threads = 1),
        # expect_error(igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), fasta = "", threads = 1),
        "The fasta file directory NA does not exist."
    )
})

test_that("returns a data.frame object", {
    result <- igblast(database = system.file("extdata/test_fasta/KIMDB_rm", package = "scifer"), system.file("extdata/test_fasta/test_igblast.txt", package = "scifer"), threads = 1)
    # expect_true(is.data.frame(result))
    expect_s3_class(result, "data.frame")
})
