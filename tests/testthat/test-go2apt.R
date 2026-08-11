
res <- go2apt("GO:2001302") # Output is AptNames, not SeqIds

test_that("`go2apt()` returns the expected result with default args", {
    expect_equal(res, c("seq.12422.143", "seq.13543.7", "seq.19617.5"))
    expect_all_true(SomaDataIO::is.AptName(res))
    expect_message(
        go2apt("GO:2001302", verbose = TRUE), 
        "`col_meta_df` not provided, using 11K annotations by default."
        )
})

test_that("`go2apt()` errors when a GO ID is not recognized", {
    expect_error(go2apt("GO:1234567"), "The provided GO term 'GO:1234567' was not found.")
    expect_error(go2apt("Pathway name string"), "The provided GO term 'Pathway name string' was not found.")
})

test_that("`go2apt()` produces expected results when `col_meta_df` is provided", {
    anno_5k <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
    res_5k <- go2apt("GO:0034975", col_meta_df = anno_5k)
    
    # Order of returned apts doesn't matter
    expect_setequal(res_5k, c("seq.5264.65", "seq.8834.58", "seq.4719.58", 
                              "seq.16588.10", "seq.4278.14", "seq.6393.63", 
                              "seq.8297.8", "seq.4959.2"))
    
    # GO:0034975 contains AptNames that are not present in the 5K menu.
    # 11450-110 should be present in 11K results, but not 5K
    expect_true("seq.11450.110" %in% setdiff(go2apt("GO:0034975"), res_5k))
})

test_that("`go2apt()` returns unique AptNames (no duplicates)", {
    res <- go2apt("GO:2001302")
    expect_false(any(duplicated(res)))
})

test_that("`go2apt()` drops AptNames with no match in col_meta_df", {
    # Provide metadata for only the 5K panel — some 11K AptNames for the GO
    # term won't be in the 5K, so the result should be a subset of 11K result
    anno_5k  <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
    res_5k   <- go2apt("GO:0034975", col_meta_df = anno_5k)
    res_11k   <- go2apt("GO:0034975")
    
    expect_true(length(res_5k) < length(res_11k))
    expect_true(all(res_5k %in% res_11k))
})

test_that("`go2apt()` can be applied over a vector of GO terms via sapply", {
    go_terms <- c("GO:2001256", "GO:0000462", "GO:2001259")
    res_list <- suppressMessages(
        sapply(go_terms, go2apt, USE.NAMES = TRUE)
    )
    
    expect_type(res_list, "list")
    expect_named(res_list, go_terms)
    expect_true(all(sapply(res_list, function(x) all(SomaDataIO::is.AptName(x)))))
    # Each term should return at least one AptName
    expect_true(all(lengths(res_list) > 0L))
})
