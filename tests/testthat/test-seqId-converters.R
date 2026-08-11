# ------ Setup

# Simple mock metadata for tests, with instances of multimapping
mock_meta <- data.frame(
    AptName = c("seq.1111.11", "seq.2222.22", "seq.3333.33", "seq.4444.44", 
                "seq.5555.55", "seq.6666.66", "seq.7777.77", "seq.8888.88"),
    EntrezGeneSymbol = c("Gene1", "Gene2", "Gene3", "Gene4 Gene5", 
                         "Gene6|Gene7|Gene8", "Gene9", "Gene9", "Gene9"),
    EntrezGeneID = c("100", "200", "300", "400 500", 
                     "600|700|800", "900", "900", "900")
)

# Real-world metadata for default argument behavior tests
real_meta <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
test_apts <- withr::with_seed(123, sample(real_meta$AptName, 3L))


# ------ Testing

#################
# Apt -> Gene
#################

test_that("`apt2gene()` uses 11K ADAT annotations by default", {
    # Verbosity is interactive-dependent, must be set to TRUE to work here
    expect_message(apt2gene(test_apts, verbose = TRUE), "`col_meta_df` not provided, using 11K annotations by default.")
    expect_equal(apt2gene(test_apts), c("SRGN", "CHEK2", "SH3GL3"))
})

test_that("`apt2gene()` converts single AptName to gene", {
    res <- apt2gene("seq.1111.11", col_meta_df = mock_meta)
    expect_equal(res, "Gene1")
})

test_that("`apt2gene()` converts a vector of AptNames to genes", {
    res <- apt2gene(c("seq.1111.11", "seq.2222.22", "seq.3333.33"), mock_meta)
    expect_equal(res, c("Gene1", "Gene2", "Gene3"))
})

test_that("`apt2gene()` errors if input is not in AptName format", {
    expect_error(apt2gene("12345-67", mock_meta),
                 "At least some values in `x` are not 'AptNames': ")
    expect_error(apt2gene("Gene1", mock_meta),
                 "At least some values in `x` are not 'AptNames': ")
})

test_that("`apt2gene()` errors on blank or NA values", {
    expect_error(apt2gene(c("seq.1111.11", ""), mock_meta),
                 "At least some values in `x` are not 'AptNames': ")
    expect_error(apt2gene(c("seq.1111.11", NA_character_), mock_meta),
                 "At least some values in `x` are not 'AptNames': NA")
})

test_that("`apt2gene()` ignores NULL values in input vector", {
    expect_equal(apt2gene(c("seq.1111.11", "seq.2222.22", NULL), mock_meta),
                 c("Gene1", "Gene2"))
})

test_that("`apt2gene()` returns NA when AptName not found in metadata", {
    expect_equal(apt2gene("seq.9999.99", mock_meta), NA_character_)
    expect_equal(apt2gene(c("seq.1111.11", "seq.9999.99"), mock_meta),
                 c("Gene1", NA_character_))
})

test_that("`apt2gene()` preserves duplicates in input", {
    res <- apt2gene(c("seq.1111.11", "seq.2222.22", "seq.1111.11"), mock_meta)
    expect_equal(res, c("Gene1", "Gene2", "Gene1"))
    expect_length(res, 3L)
})

test_that("`apt2gene(collapse = TRUE)` keeps multi-gene mappings as single string", {
    # Space-delimited genes
    expect_equal(apt2gene("seq.4444.44", mock_meta, collapse = TRUE), "Gene4 Gene5")
    # Pipe-delimited genes
    expect_equal(apt2gene("seq.5555.55", mock_meta, collapse = TRUE), "Gene6|Gene7|Gene8")
})

test_that("`apt2gene(collapse = FALSE)` splits multi-gene mappings", {
    expect_equal(apt2gene("seq.4444.44", mock_meta, collapse = FALSE),
                 c("Gene4", "Gene5"))
    expect_equal(apt2gene("seq.5555.55", mock_meta, collapse = FALSE),
                 c("Gene6", "Gene7", "Gene8"))
})

test_that("`apt2gene()` uses different gene ID type when specified", {
    res <- apt2gene("seq.1111.11", mock_meta, id_type = "EntrezGeneID")
    expect_equal(res, "100")
})

test_that("`apt2gene()` errors with list-specific message when a list is provided", {
    # List check must come before is.AptName() check, otherwise wrong error fires
    expect_error(apt2gene(list("seq.1111.11", "seq.2222.22"), mock_meta),
                 "`x` must be a vector, not a list")
})

test_that("`apt2gene()` errors when `col_meta_df` is not a data frame", {
    expect_error(
        apt2gene("seq.1111.11", col_meta_df = "not_a_df"),
        "`col_meta_df` must be a data frame, tibble, or similar."
    )
    expect_error(
        apt2gene("seq.1111.11", col_meta_df = list(AptName = "seq.1111.11")),
        "`col_meta_df` must be a data frame, tibble, or similar."
    )
})


#################
# Gene -> Apt
#################

test_that("`gene2apt()` converts single gene to AptName", {
    res <- gene2apt("Gene1", mock_meta)
    expect_equal(res, "seq.1111.11")
})

test_that("`gene2apt()` converts a vector of genes to AptNames", {
    res_vec <- gene2apt(c("Gene1", "Gene2", "Gene3"), mock_meta)
    expect_equal(res_vec, c("seq.1111.11", "seq.2222.22", "seq.3333.33"))
})

test_that("`gene2apt()` uses specified 'id_type'", {
    res_id <- gene2apt("100", mock_meta, id_type = "EntrezGeneID")
    res_id_na <- gene2apt("GENE1", mock_meta, id_type = "EntrezGeneID")
    expect_equal(res_id, "seq.1111.11")
    expect_true(is.na(res_id_na))
})

test_that("`gene2apt()` drops blank, NULL, and NA values in input", {
    expect_equal(gene2apt(c("Gene1", "Gene2", NA), mock_meta),
                 c("seq.1111.11", "seq.2222.22"))
    expect_equal(gene2apt(c("Gene1", "Gene2", ""), mock_meta),
                 c("seq.1111.11", "seq.2222.22"))
    expect_equal(gene2apt(c("Gene1", "Gene2", NULL), mock_meta),
                 c("seq.1111.11", "seq.2222.22"))
})

test_that("`gene2apt()` returns NA when gene not found in metadata", {
    na_res <- gene2apt("NA_TEST", mock_meta)
    na_vec <- gene2apt(c("Gene1", "NA_TEST"), mock_meta)
    expect_equal(na_res, NA_character_)
    expect_equal(na_vec, c("seq.1111.11", NA_character_))
})

test_that("`gene2apt()` preserves duplicates in input", {
    dup_res <- gene2apt(c("Gene1", "Gene2", "Gene1"), mock_meta)
    expect_equal(dup_res, c("seq.1111.11", "seq.2222.22", "seq.1111.11"))
    expect_length(dup_res, 3L)
})

test_that("`gene2apt(collapse = TRUE)` collapses multiple AptNames for one gene", {
    result <- gene2apt("Gene9", mock_meta, collapse = TRUE)
    expect_equal(result, "seq.6666.66|seq.7777.77|seq.8888.88")
})

test_that("`gene2apt(collapse = TRUE)` uses custom separator", {
    result <- gene2apt("Gene9", mock_meta, collapse = TRUE, sep = ", ")
    expect_equal(result, "seq.6666.66, seq.7777.77, seq.8888.88")
})

test_that("`gene2apt(collapse = FALSE)` returns all multimapped AptNames individually", {
    result <- gene2apt("Gene9", mock_meta, collapse = FALSE)
    expect_equal(result, c("seq.6666.66", "seq.7777.77", "seq.8888.88"))
})

test_that("`gene2apt(collapse = FALSE)` handles mixed single/multi mappings", {
    result <- gene2apt(c("Gene1", "Gene9"), mock_meta, collapse = FALSE)
    expect_equal(result, c("seq.1111.11", "seq.6666.66", "seq.7777.77", "seq.8888.88"))
})

test_that("`gene2apt()` errors if list is provided as input", {
    expect_error(gene2apt(list("Gene1", "Gene2"), mock_meta),
                 "`x` must be a vector, not a list")
})


# ----- Fuzzy column name resolution

test_that("`apt2gene()` works when col_meta_df has dot-separated column name", {
    mock_meta_ipp <- mock_meta
    names(mock_meta_ipp)[names(mock_meta_ipp) == "EntrezGeneSymbol"] <- "Entrez.Gene.Symbol"
    res_ipp   <- apt2gene("seq.1111.11", col_meta_df = mock_meta_ipp)
    res_array <- apt2gene("seq.1111.11", col_meta_df = mock_meta)
    expect_equal(res_ipp, res_array)
})

test_that("`gene2apt()` works when col_meta_df has dot-separated column name", {
    mock_meta_ipp <- mock_meta
    names(mock_meta_ipp)[names(mock_meta_ipp) == "EntrezGeneSymbol"] <- "Entrez.Gene.Symbol"
    res_ipp   <- gene2apt("Gene1", col_meta_df = mock_meta_ipp)
    res_array <- gene2apt("Gene1", col_meta_df = mock_meta)
    expect_equal(res_ipp, res_array)
})

