test_that("output dataframe 96 well", {
    suppressMessages(
        index_sort_data <- fcs_processing(
            folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
            compensation = TRUE, plate_wells = 96,
            probe1 = "Pre.F", probe2 = "Post.F",
            posvalue_probe1 = 600, posvalue_probe2 = 400
        )
    )
    expect_type(index_sort_data, "list")
    expect_s3_class(index_sort_data[[1]], "data.frame")
    expect_s3_class(index_sort_data[[2]], "data.frame")
})

test_that("output dataframe 384-well", {
    suppressMessages(
        index_sort_data <- fcs_processing(
            folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
            compensation = TRUE, plate_wells = 384,
            probe1 = "Pre.F", probe2 = "Post.F",
            posvalue_probe1 = 600, posvalue_probe2 = 400
        )
    )
    expect_type(index_sort_data, "list")
    expect_s3_class(index_sort_data[[1]], "data.frame")
    expect_s3_class(index_sort_data[[2]], "data.frame")
})

test_that("applying compensation", {
    suppressMessages(
        expect_message(
            fcs_processing(
                folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
                compensation = TRUE, plate_wells = 384,
                probe1 = "Pre.F", probe2 = "Post.F",
                posvalue_probe1 = 600, posvalue_probe2 = 400
            ),
            "Samples were compensated using the compensation saved on fsc index files."
        )
    )

    suppressMessages(
        expect_message(
            fcs_processing(
                folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
                compensation = FALSE, plate_wells = 384,
                probe1 = "Pre.F", probe2 = "Post.F",
                posvalue_probe1 = 600, posvalue_probe2 = 400
            ),
            "Samples were not compensated."
        )
    )
    expect_error(
        fcs_processing(
            folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
            compensation = 124, plate_wells = 96
        ),
        "Compensation argument should be TRUE or FALSE."
    )
})

test_that("only allow plates with 96 or 384-wells", {
    suppressMessages(
        expect_error(
            fcs_processing(
                folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
                compensation = TRUE, plate_wells = 123
            ),
            "Only 96 or 384-well plates are supported"
        )
    )
    suppressMessages(
        expect_error(
            fcs_processing(
                folder_path = system.file("/extdata/fcs_index_sorting", package = "scifer"),
                compensation = TRUE, plate_wells = "string"
            ),
            "Only 96 or 384-well plates are supported"
        )
    )
})
