
res <- path2apt("GO:2001302", verbose = FALSE)

test_that("`path2apt()` returns the expected result with default args", {
    expect_type(res, "list")
    expect_named(res, "GO:2001302")
    expect_equal(res[[1L]], c("seq.12422.143", "seq.13543.7", "seq.19617.5"))
    expect_all_true(SomaDataIO::is.AptName(res[[1L]]))
})

test_that("`path2apt()` errors when `col_meta_df` is not a data frame", {
    list_meta <- list(t_tests$AptName, t_tests$EntrezGeneSymbol)
    
    expect_error(
        path2apt("GO:2001302", col_meta_df = list_meta), 
        "`col_meta_df` must be a data frame or tibble.", 
        fixed = TRUE
    )
})

test_that("`path2apt()` utilizes user-provided `col_meta_df` as expected", {
    path <- pathway_map[pathway_map$pathway_id == "GO:2001302", ]
    
    # Create dummy gene symbol -> AptName mappings
    test_meta <- data.frame(EntrezGeneSymbol = c("ALOX12", "CYP4F3", "PTGR1"),
                            AptName = c("seq.test1.01", "seq.test2.02", "seq.test3.03"))
    res_list <- list(`GO:2001302` = c("seq.test1.01", "seq.test2.02", "seq.test3.03"))
    
    # Should return my 3 dummy AptNames
    expect_equal(res_list, path2apt("GO:2001302", col_meta_df = test_meta))
})

test_that("`path2apt()` utilizes default `col_meta_df` as expected", {
    col_meta <- getAnalyteInfo(example_data_11k)
    genes <- dplyr::filter(pathway_map, pathway_id == "GO:2001302")
    map <- dplyr::filter(col_meta, EntrezGeneSymbol %in% genes$gene_symbol)$AptName
    res_list <- list(`GO:2001302` = map)
    
    expect_message(path2apt("GO:2001302", verbose = TRUE),
                   "`col_meta_df` not provided, using 11K annotations by default.",
                   fixed = TRUE)
    expect_equal(res_list, path2apt("GO:2001302", col_meta_df = col_meta))
})

test_that("`path2apt()` errors when a pathway ID is not recognized", {
    expect_error(path2apt("GO:1234567"),
                 "Pathway ID(s) not found in `pathway_map`: GO:1234567",
                 fixed = TRUE)
    expect_error(path2apt(c("GO:2001302", "Pathway name string", "M0000")),
                 "Pathway ID(s) not found in `pathway_map`: Pathway name string, M0000",
                 fixed = TRUE)
})

test_that("`path2apt()` errors when `x` is not a character vector", {
    expect_error(path2apt(list("GO:2001302")), "`x` must be a character vector.")
    expect_error(path2apt(123), "`x` must be a character vector.")
})

test_that("`path2apt()` matches `go2apt()` when GO terms are provided as input", {
    go_terms <- c("GO:2001256", "GO:0000462", "GO:2001259")
    res_list <- path2apt(go_terms, verbose = FALSE)
    
    lapply(go_terms, function(x) {
        expect_equal(res_list[[x]], go2apt(x, verbose = FALSE))
    })
})

test_that("`path2apt()` works with non-GO pathways", {
    h_genes <- pathway_map[pathway_map$pathway_id == "M5890", ]$gene_symbol
    res_h   <- path2apt("M5890", verbose = FALSE)
    
    # Order of returned apts doesn't matter
    expect_setequal(res_h[["M5890"]], gene2apt(h_genes, collapse = FALSE))
})

test_that("`path2apt()` produces expected results when `col_meta_df` is provided", {
    anno_5k <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
    res_5k <- path2apt("GO:0034975", col_meta_df = anno_5k)
    
    # Order of returned apts doesn't matter
    expect_setequal(res_5k[[1L]], c("seq.5264.65", "seq.8834.58", "seq.4719.58", 
                                    "seq.16588.10", "seq.4278.14", "seq.6393.63", 
                                    "seq.8297.8", "seq.4959.2"))
    
    # GO:0034975 contains AptNames that are not present in the 5K menu.
    # 11450-110 should be present in 11K results, but not 5K
    res_11k <- path2apt("GO:0034975", verbose = FALSE)
    expect_true("seq.11450.110" %in% setdiff(res_11k[[1L]], res_5k[[1L]]))
})

test_that("`path2apt()` works with `id_type = 'EntrezGeneID'`", {
    h_ids <- pathway_map[pathway_map$pathway_id == "M5890", ]$entrez_id
    res_e <- path2apt("M5890", id_type = "EntrezGeneID", verbose = FALSE)
    
    expect_all_true(SomaDataIO::is.AptName(res_e[[1L]]))
    expect_setequal(res_e[["M5890"]], 
                    gene2apt(h_ids, id_type = "EntrezGeneID", collapse = FALSE))
})

test_that("`path2apt()` returns unique AptNames (no duplicates)", {
    res_list <- path2apt(c("GO:2001302", "M5890"), verbose = FALSE)
    expect_false(any(vapply(res_list, anyDuplicated, integer(1)) > 0L))
})

test_that("`path2apt()` preserves input order and drops duplicate IDs", {
    ids      <- c("M5890", "GO:2001256", "GO:0000462", "M5890")
    res_list <- path2apt(ids, verbose = FALSE)
    
    expect_named(res_list, unique(ids))
    expect_length(res_list, 3L)
})

test_that("`path2apt()` returns `character(0)` for pathways with no matching analytes", {
    anno_min <- data.frame(AptName = "seq.1234.56", EntrezGeneSymbol = "NOTAGENE")
    res_none <- path2apt(c("GO:2001302", "M5890"), col_meta_df = anno_min)
    
    expect_named(res_none, c("GO:2001302", "M5890"))
    expect_equal(lengths(res_none), c("GO:2001302" = 0L, "M5890" = 0L))
    expect_type(res_none[[1L]], "character")
})
