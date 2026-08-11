# Setup ------

# Mock column metadata for unit tests
mock_meta <- data.frame(
    AptName = c("seq.1111.11", "seq.2222.22", "seq.3333.33", "seq.4444.44",
                "seq.5555.55", "seq.6666.66", "seq.7777.77", "seq.8888.88"),
    EntrezGeneSymbol = c("Gene1", "Gene2", "Gene3", "Gene4 Gene5",
                         "Gene6|Gene7|Gene8", "Gene9", "Gene9", "Gene9"),
    EntrezGeneID = c("100", "200", "300", "400 500",
                     "600|700|800", "900", "900", "900")
)

# Mock rankings (values = simulated fold change)
# Gene9 has 3 analytes: seq.6666.66, seq.7777.77, seq.8888.88
mock_ranks <- c(seq.6666.66 = 5.0, seq.1111.11 = 4.5, seq.7777.77 = 4.0,
                seq.2222.22 = 3.5, seq.3333.33 = 3.0, seq.8888.88 = 2.5,
                seq.4444.44 = 2.0, seq.5555.55 = 1.5)

# Mock pathways
mock_path_simple   <- list(Simple_Pathway = c("Gene1", "Gene2", "Gene3"))
mock_path_multimer <- list(Multimer_Pathway = c("Gene1", "Gene9"))  # Gene9 has 3 aptamers
mock_path_missing  <- list(Missing_Pathway = c("Gene1", "Test1", "Gene2"))
mock_path_id       <- list(EntrezID_Pathway = c("100", "200", "900"))
mock_path_empty    <- list(Empty_Pathway = c("Test1", "Test2"))
mock_path_multiple <- list(Pathway1 = c("Gene1", "Gene2"),
                           Pathway2 = c("Gene3", "Gene9"),
                           Pathway3 = c("Gene1", "Gene9", "Gene2"))

# Real-world data 
ranks <- abs(t_tests$log2_fc)
names(ranks) <- t_tests$AptName
ranks <- sort(ranks, decreasing = TRUE)
meta <- getAnalyteInfo(example_data_11k)


# Testing -------
test_that("`genePath2aptPath()` returns expected list with defaults", {
    out <- genePath2aptPath(path_list = mock_path_simple,
                            col_meta_df = mock_meta)

    expect_type(out, "list")
    expect_length(out, 1L)
    expect_named(out, "Simple_Pathway")
})

test_that("`genePath2aptPath()` converts genes to AptNames correctly", {
    out <- genePath2aptPath(path_list = mock_path_simple,
                            col_meta_df = mock_meta)

    expect_equal(sort(unname(out$Simple)),
                 c("seq.1111.11", "seq.2222.22", "seq.3333.33"))
    expect_equal(sort(names(out$Simple)), c("Gene1", "Gene2", "Gene3"))
})

test_that("`genePath2aptPath()` retains all aptamers when gene maps to multiple", {
    out <- genePath2aptPath(path_list = mock_path_multimer,
                            col_meta_df = mock_meta)

    # Gene1 (1 apt) + Gene9 (3 apts) = 4 total
    expect_length(out$Multimer_Pathway, 4L)
    # All 3 Gene9 aptamers should be present
    gene9_apts <- unname(out$Multimer_Pathway[names(out$Multimer_Pathway) == "Gene9"])
    expect_setequal(gene9_apts, c("seq.6666.66", "seq.7777.77", "seq.8888.88"))
    expect_equal(out$Multimer_Pathway["Gene1"], c(Gene1 = "seq.1111.11"))
})

test_that("`genePath2aptPath()` drops genes not found in metadata", {
    out <- genePath2aptPath(path_list = mock_path_missing,
                            col_meta_df = mock_meta)

    expect_length(out$Missing, 2L)
    expect_true("Gene1" %in% names(out$Missing))
    expect_true("Gene2" %in% names(out$Missing))
    expect_false("GeneNotFound" %in% names(out$Missing))
})

test_that("`genePath2aptPath()` returns empty named vector when no genes match", {
    out <- genePath2aptPath(path_list = mock_path_empty,
                            col_meta_df = mock_meta)

    expect_length(out$Empty, 0L)
    expect_type(out$Empty, "character")
})

test_that("`genePath2aptPath()` handles EntrezGeneID id_type", {
    out <- genePath2aptPath(path_list = mock_path_id,
                            col_meta_df = mock_meta,
                            id_type = "EntrezGeneID")

    # ID 100 (1 apt) + 200 (1 apt) + 900 (3 apts) = 5 total
    expect_length(out$EntrezID_Pathway, 5L)
    # All 3 aptamers for ID 900 (Gene9) should be present
    id900_apts <- unname(out$EntrezID_Pathway[names(out$EntrezID_Pathway) == "900"])
    expect_setequal(id900_apts, c("seq.6666.66", "seq.7777.77", "seq.8888.88"))
})

test_that("`genePath2aptPath()` processes >1 pathways correctly", {
    out <- genePath2aptPath(path_list = mock_path_multiple,
                            col_meta_df = mock_meta,
                            verbose = FALSE)

    expect_length(out, 3L)
    expect_named(out, c("Pathway1", "Pathway2", "Pathway3"))
    expect_length(out$Pathway1, 2L)  # Gene1 + Gene2
    expect_length(out$Pathway2, 4L)  # Gene3 + Gene9 (3 apts)
    expect_length(out$Pathway3, 5L)  # Gene1 + Gene9 (3 apts) + Gene2
})

test_that("`genePath2aptPath()` preserves pathway independence", {
    # Each pathway should be processed independently
    out <- genePath2aptPath(path_list = mock_path_multiple,
                            col_meta_df = mock_meta,
                            verbose = FALSE)

    # Gene1 appears in Path1 and Path3; both should have the same result
    expect_equal(out$Pathway1["Gene1"], out$Pathway3["Gene1"])
    # Gene9 appears in Path2 and Path3; all 3 aptamers should be present in both
    gene9_p2 <- sort(unname(out$Pathway2[names(out$Pathway2) == "Gene9"]))
    gene9_p3 <- sort(unname(out$Pathway3[names(out$Pathway3) == "Gene9"]))
    expect_equal(gene9_p2, gene9_p3)
    expect_setequal(gene9_p2, c("seq.6666.66", "seq.7777.77", "seq.8888.88"))
})

test_that("`genePath2aptPath()` handles genes from multi-gene aptamer entries", {
    # Gene4 and Gene5 share the same aptamer (seq.4444.44)
    path_shared <- list(Shared = c("Gene4", "Gene5"))
    out <- genePath2aptPath(path_list = path_shared,
                            col_meta_df = mock_meta)

    expect_length(out$Shared, 2L)
    # Both genes should map to the same aptamer
    expect_equal(unname(out$Shared["Gene4"]), "seq.4444.44")
    expect_equal(unname(out$Shared["Gene5"]), "seq.4444.44")
})

test_that("`genePath2aptPath()` deduplicates genes within a pathway", {
    path_dups <- list(Dups = c("Gene1", "Gene2", "Gene1", "Gene1"))
    out <- genePath2aptPath(path_list = path_dups,
                            col_meta_df = mock_meta)

    # Duplicates should be removed, resulting in 2 unique genes
    expect_length(out$Dups, 2L)
    expect_equal(sort(names(out$Dups)), c("Gene1", "Gene2"))
})

test_that("`genePath2aptPath()` works when col_meta_df has dot-separated column name", {
    mock_meta_ipp <- mock_meta
    names(mock_meta_ipp)[names(mock_meta_ipp) == "EntrezGeneSymbol"] <- "Entrez.Gene.Symbol"
    res_ipp <- genePath2aptPath(path_list = mock_path_simple,
                                col_meta_df = mock_meta_ipp,
                                verbose = FALSE)
    res_array <- genePath2aptPath(path_list = mock_path_simple,
                                  col_meta_df = mock_meta,
                                  verbose = FALSE)
    expect_equal(res_ipp, res_array)
})



