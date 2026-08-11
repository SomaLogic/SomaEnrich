# Setup ----
meta      <- SomaDataIO::getAnalyteInfo(example_data_11k)
gmt_file  <- system.file("extdata", "msigdb_h_androgen_response.gmt",
                         package = "SomaEnrich")
gmt_file2 <- system.file("extdata", "msigdb_c2_reactome.gmt",
                         package = "SomaEnrich")


# Testing ----
test_that("`gmt_gene2apt()` returns a named list with default args", {
    res <- gmt_gene2apt(gmt_gene_file = gmt_file,
                        col_meta_df   = meta)

    expect_type(res, "list")
    expect_named(res)
    expect_gt(length(res), 0L)
})

test_that("`gmt_gene2apt()` output contains AptNames", {
    res <- gmt_gene2apt(gmt_gene_file = gmt_file,
                        col_meta_df   = meta)

    # Every non-empty pathway should contain valid AptNames
    non_empty <- res[lengths(res) > 0L]
    expect_true(all(SomaDataIO::is.AptName(unlist(non_empty))))
})

test_that("`gmt_gene2apt()` preserves pathway names from the GMT file", {
    original_paths <- names(fgsea::gmtPathways(gmt_file))
    res <- gmt_gene2apt(gmt_gene_file = gmt_file,
                        col_meta_df   = meta)

    expect_true(all(names(res) %in% original_paths))
})

test_that("`gmt_gene2apt()` respects `id_type = 'EntrezGeneID'`", {
    res_sym <- gmt_gene2apt(gmt_gene_file = gmt_file,
                            col_meta_df   = meta,
                            id_type       = "EntrezGeneSymbol")
    res_id  <- gmt_gene2apt(gmt_gene_file = gmt_file,
                            col_meta_df   = meta,
                            id_type       = "EntrezGeneID")

    # Both should return valid AptNames; results may differ slightly due to
    # ID coverage differences but structure should be identical
    expect_type(res_id, "list")
    expect_named(res_id)
    expect_equal(names(res_sym), names(res_id))
})

test_that("`gmt_gene2apt()` result matches manual `genePath2aptPath()` call", {
    gene_paths <- fgsea::gmtPathways(gmt_file)
    manual_res <- genePath2aptPath(gene_paths, col_meta_df = meta)
    gmt_res    <- gmt_gene2apt(gmt_gene_file = gmt_file, col_meta_df = meta)

    expect_equal(names(gmt_res), names(manual_res))
    # Values should match (names within each path vector may differ in order)
    expect_true(all(mapply(function(a, b) setequal(unname(a), unname(b)),
                           gmt_res, manual_res)))
})

test_that("`gmt_gene2apt(write = TRUE)` writes a file when path is provided", {
    out_file <- tempfile(fileext = ".gmt")
    on.exit(unlink(out_file), add = TRUE)

    skip_on_os("windows")
    
    # This message doesn't print as expected on Windows, but will pass on macOS
    expect_message(
        gmt_gene2apt(gmt_gene_file = gmt_file,
                     col_meta_df   = meta,
                     write         = TRUE,
                     gmt_apt_file  = out_file,
                     verbose       = TRUE),
        paste("Writing converted GMT file to", out_file)
    )
    expect_true(file.exists(out_file))
    expect_gt(file.size(out_file), 0L)
})

test_that("`gmt_gene2apt(write = TRUE)` errors when no path is provided", {
  expect_error(
    gmt_gene2apt(gmt_gene_file = gmt_file,
                 col_meta_df   = meta,
                 write         = TRUE,
                 verbose       = FALSE), 
    "Please provide a filename to `gmt_apt_file` when `write = TRUE`.")
})

test_that("`gmt_gene2apt()` works on a multi-pathway GMT file", {
    res <- gmt_gene2apt(gmt_gene_file = gmt_file2,
                        col_meta_df   = meta)

    expect_type(res, "list")
    expect_gt(length(res), 1L)  # More than one pathway
})

