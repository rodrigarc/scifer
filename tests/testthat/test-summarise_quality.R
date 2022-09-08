test_that("Output list with df and list", {
    suppressMessages(
        sf <- summarise_quality(
            folder_sequences = system.file("extdata/sorted_sangerseq", package = "scifer"),
            secondary.peak.ratio = 0.33,
            trim.cutoff = 0.01,
            processor = 1
        )
    )
    expect_type(sf, "list")
    expect_type(sf[["quality_scores"]], "list")
    expect_s3_class(sf[["summaries"]], "data.frame")
})

test_that("Input", {
    suppressMessages(
        expect_error(summarise_quality())
    )
    suppressMessages(
        expect_error(summarise_quality("non_exitent_path"))
    )
})
