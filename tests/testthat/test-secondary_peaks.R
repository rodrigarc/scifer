test_that("Output list with df and s4", {
    s4_sangerseq <- sangerseqR::readsangerseq(
        system.file("/extdata/sorted_sangerseq/E18_C1/A1_3_IgG_Inner.ab1", package = "scifer")
    )
    processed_seq <- secondary_peaks(s4_sangerseq)

    expect_type(processed_seq, "list")
    expect_s3_class(processed_seq[["secondary.peaks"]], "data.frame")
    expect_s4_class(processed_seq[["read"]], "sangerseq")
})
