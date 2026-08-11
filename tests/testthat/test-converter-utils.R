#----- Setup

meta <- SomaDataIO::getAnalyteInfo(example_data_11k)
genes <- withr::with_seed(303, sample(meta$EntrezGeneSymbol, 3))
genes_mm <- c("KLK3", "CGA")
apts <- withr::with_seed(303, sample(meta$AptName, 3)) # Must use same seed as genes!
ranks_unsrt <- c(GENE_A = 0.01, GENE_A = 0.05, GENE_A = 0.09, GENE_B = 0.02, GENE_C = 0.03)
ranks_srt <- sort(ranks_unsrt)

# Total unique pairs of SOMAmer reagent + genes in the 7K menu
pairs_7k <- meta |>
    dplyr::select(AptName, EntrezGeneSymbol) |>
    tidyr::separate_rows(EntrezGeneSymbol, sep = "\\,|\\s+|\\|") |>
    dplyr::filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "") |>
    dplyr::mutate(Key = paste0(AptName, ",", EntrezGeneSymbol)) |>
    unique()


# ----- Testing

####### .make_id_map() ########
test_that(".make_id_map() returns expected object with defaults", {
    res <- .make_id_map(meta)
    res$Key <- paste0(res$AptName, ",", res$EntrezGeneSymbol)
    
    expect_true(inherits(res, "data.frame"))
    expect_equal(nrow(res), nrow(pairs_7k))
    expect_named(res, c("AptName", "EntrezGeneSymbol", "Key"))
    expect_false(any(is.na(res$EntrezGeneSymbol)))
    expect_false(any(res$EntrezGeneSymbol == ""))
    expect_length(setdiff(res$Key, pairs_7k$Key), 0L)
    expect_true(all(res$AptName %in% pairs_7k$AptName))
    expect_true(all(res$EntrezGeneSymbol %in% pairs_7k$EntrezGeneSymbol))
})

test_that(".make_id_map() splits '|'-delimited IDs, by default", {
    res <- .make_id_map(meta)
    d <- setdiff(meta$EntrezGeneSymbol, res$EntrezGeneSymbol)
    
    # AptNames aren't delimited and should therefore be identical
    expect_true(all(res$AptName %in% meta$AptName))
    # 1st character is a "" and can be removed for this test
    expect_true(all(grepl("\\|", d[-1]))) 
})

test_that(".make_id_map() filters metadata when genes are provided", {
    res <- .make_id_map(meta, x = genes)
    genes_srt <- sort(res$EntrezGeneSymbol)
    
    expect_true(inherits(res, "data.frame"))
    expect_equal(nrow(res), 4L) # 2 analytes for gene TAGLN2
    expect_equal(genes_srt, c("RABL3", "SOD2", "TAGLN2", "TAGLN2"))
})

test_that(".make_id_map() filters metadata when AptNames are provided", {
    res <- .make_id_map(meta, x = apts)
    
    expect_true(inherits(res, "data.frame"))
    expect_equal(nrow(res), 3L) # One analyte per gene
    expect_equal(sort(res$EntrezGeneSymbol), 
                 c("RABL3", "SOD2", "TAGLN2"))
})

test_that(".make_id_map() returns empty data frame when x matches nothing", {
    res <- .make_id_map(meta, x = "TEST")
    
    expect_true(inherits(res, "data.frame"))
    expect_equal(nrow(res), 0L)
    expect_named(res, c("AptName", "EntrezGeneSymbol"))
})

test_that(".make_id_map() uses different gene_col when specified", {
    res <- .make_id_map(meta, gene_col = "EntrezGeneID")
    split_ids <- unique(unlist(strsplit(meta$EntrezGeneID, "\\||\\s")))
    
    expect_named(res, c("AptName", "EntrezGeneID"))
    expect_true(all(res$EntrezGeneID %in% split_ids))
})

test_that(".make_id_map() uses custom apt_col when specified", {
    res <- .make_id_map(meta, apt_col = "SeqId")
    
    expect_named(res, c("SeqId", "EntrezGeneSymbol"))
    expect_true(all(SomaDataIO::is.SeqId(res$SeqId)))
})

test_that(".make_id_map() expands multi-gene entries into separate rows", {
    res <- .make_id_map(meta, x = genes_mm)
    
    # Each of these multimers maps to multiple individual genes after expansion
    expect_true(nrow(res) > length(genes_mm))
    expect_equal(sum(grepl("CGA", res$EntrezGeneSymbol)), 5L)
    expect_equal(sum(grepl("KLK3", res$EntrezGeneSymbol)), 4L)
    expect_true(all(genes_mm %in% res$EntrezGeneSymbol))
})

test_that(".make_id_map() removes empty and NA gene entries", {
    # Create a mock df with empty/NA gene values
    mock_meta <- data.frame(AptName = c("seq.12345.12", "seq.987654.98", "seq.321321.32", "seq.44444.44"),
                            EntrezGeneSymbol = c("GENE1", "", NA_character_, "GENE2"))
    res <- .make_id_map(mock_meta)
    
    expect_equal(nrow(res), 2L)
    expect_equal(sort(res$EntrezGeneSymbol), c("GENE1", "GENE2"))
    expect_false(any(is.na(res$EntrezGeneSymbol)))
    expect_false(any(res$EntrezGeneSymbol == ""))
})

test_that(".make_id_map() removes duplicate rows", {
    mock_meta <- data.frame(
        AptName = c("seq.1234.56", "seq.1234.56", "seq.9876.54"),
        EntrezGeneSymbol = c("GENE1", "GENE1", "GENE2")
    )
    res <- .make_id_map(mock_meta)
    
    expect_equal(nrow(res), 2L)  # Duplicate rows should be removed
})

test_that(".make_id_map() handles empty input data frame", {
    empty_meta <- data.frame(AptName = character(0), EntrezGeneSymbol = character(0))
    res <- .make_id_map(empty_meta)
    
    expect_true(inherits(res, "data.frame"))
    expect_equal(nrow(res), 0L)
})


####### .collapse_by_gene() ########
test_that(".collapse_by_gene() collapses AptNames per gene using default separator", {
    # APOE maps to multiple AptNames
    id_map <- .make_id_map(meta, x = "APOE")
    res <- .collapse_by_gene(id_map)
    apoe_apts <- c("seq.2418.55", "seq.2937.10", "seq.2938.55", "seq.5312.49")
    
    expect_named(res, "APOE")
    expect_length(res, 1L)
    expect_true(grepl("\\|", res)) # Default separator is "|"
    
    # All APOE analytes are correctly identified and collapsed
    expect_equal(length(strsplit(res, "\\|")[[1]]), nrow(id_map))
    expect_true(all(apoe_apts %in% strsplit(res, "\\|")[[1]]))
})

test_that(".collapse_by_gene() uses custom separator when specified", {
    id_map <- .make_id_map(meta, x = "APOE")
    res <- .collapse_by_gene(id_map, sep = ";")
    matches <- gregexpr(pattern = ";", text = res, fixed = TRUE)
    
    expect_true(grepl(";", res))
    expect_length(matches[[1]], 3L) # 3 ";" characters to separate 4 APOE values
    expect_false(grepl("\\|", res))
})

test_that(".collapse_by_gene() groups values by shared gene only", {
    # Create mock with multiple genes, each having multiple aptamers
    mock_map <- data.frame(
        AptName = c("seq.1.1", "seq.2.2", "seq.3.3", "seq.4.4", "seq.5.5"),
        EntrezGeneSymbol = c("GENE_A", "GENE_A", "GENE_B", "GENE_B", "GENE_B")
    )
    res <- .collapse_by_gene(mock_map)
    
    expect_length(res, 2L)
    expect_equal(res[1], c(GENE_A = "seq.1.1|seq.2.2"))
    expect_equal(res[2], c(GENE_B = "seq.3.3|seq.4.4|seq.5.5"))
})

test_that(".collapse_by_gene() removes duplicate AptNames within a gene", {
    mock_map <- data.frame(
        AptName = c("seq.1234.56", "seq.1234.56", "seq.9876.54"), # Duplicate seq.1234.56
        EntrezGeneSymbol = c("GENE_A", "GENE_A", "GENE_A")
    )
    res <- .collapse_by_gene(mock_map)
    
    expect_length(res, 1L)
    expect_equal(res, c(GENE_A = "seq.1234.56|seq.9876.54")) # Only 2 unique AptNames
})

test_that(".collapse_by_gene() uses custom column names when specified", {
    mock_map <- data.frame(SeqId = c("12345-1", "67890-2"),
                           GeneID = c("100", "100"))
    res <- .collapse_by_gene(mock_map, apt_col = "SeqId", gene_col = "GeneID")
    
    expect_length(res, 1L)
    expect_equal(res, c(`100` = "12345-1|67890-2"))
})

test_that(".collapse_by_gene() returns NA for genes with no valid aptamers", {
    mock_map <- data.frame(
        AptName = c(NA_character_, NA_character_),
        EntrezGeneSymbol = c("GENE_X", "GENE_X")
    )
    res <- .collapse_by_gene(mock_map)
    
    expect_length(res, 1L)
    expect_true(is.na(res))
})


####### .select_best() ########
test_that(".select_best() returns expected list structure", {
    res <- .select_best(ranks_srt)
    
    expect_type(res, "list")
    expect_named(res, c("ranks", "id_list", "keep_list")) 
})

test_that(".select_best() output contains indices of input ranks", {
    res <- .select_best(ranks_srt)
    a_idx <- grep("GENE_A", names(res$ranks)) # indices in lists come from this vector
    b_idx <- grep("GENE_B", names(res$ranks))
    c_idx <- grep("GENE_C", names(res$ranks))
    
    # All indices in id_list & keep_list should match locations in this vector
    expect_true(res$keep_list$GENE_A %in% a_idx)
    expect_true(res$keep_list$GENE_B %in% b_idx)
    expect_true(res$keep_list$GENE_C %in% c_idx)
})

test_that(".select_best() tie-breaking uses first occurrence when values are equal", {
    tied_ranks <- c(GENE_B = 0.01, GENE_A = 0.05, GENE_A = 0.05)
    
    res_min <- .select_best(tied_ranks, resolve_method = "min")
        res_max <- .select_best(tied_ranks, resolve_method = "max")

    # Both GENE_A values are equal, which.min returns first match
    expect_equal(res_min$keep_list$GENE_A, 2L)
    expect_equal(res_max$keep_list$GENE_A, 2L)
})

test_that(".select_best() correctly identifies all duplicate IDs in 'ranks'", {
    ranks <- c(A = 1, B = 2, A = 3, C = 4, B = 5, A = 6)
    res <- .select_best(ranks, resolve_method = "rank")
    
    expect_named(res$id_list, c("A", "B", "C"))
    expect_equal(res$id_list$A, c(1L, 3L, 6L))
    expect_equal(res$id_list$B, c(2L, 5L))
    expect_equal(res$id_list$C, 4L)
    
    expect_equal(res$keep_list$A, 1L)
    expect_equal(res$keep_list$B, 2L)
    expect_equal(res$keep_list$C, 4L)
})

test_that(".select_best() selects largest absolute value by default", {
    # GENE_A has values -5 (pos 1) and 2 (pos 3)
    # resolve_method='abs' should select position 1 (highest absolute value = 5)
    ranks_test <- c(GENE_A = -5, GENE_B = 1, GENE_A = 2)
    res <- .select_best(ranks_test)
    
    expect_equal(res$keep_list$GENE_A, 1L)
    expect_equal(res$keep_list$GENE_B, 2L) # Sanity check 
    
    # Another example where largest negative value is retained
    ranks_neg <- c(GENE_X = 3, GENE_X = -10, GENE_Y = 5, GENE_X = 8)
    res_neg <- .select_best(ranks_neg)
    
    expect_equal(res_neg$keep_list$GENE_X, 2L) # Position 2 has highest abs value (|-10| = 10)
    
    # Example where positive value retained
    ranks_pos <- c(GENE_Z = -5, GENE_Z = 8, GENE_W = 1)
    res_pos <- .select_best(ranks_pos)
    
    expect_equal(res_pos$keep_list$GENE_Z, 2L) # Position 2 has highest abs value
})

# Ensures that pre-sorting of ranks actually affects the selection process
test_that(".select_best() keeps first occurrence when `resolve_method = 'rank'`", {
    # If sorted by ascending p-value, 1st occurrence of GENE_A should be kept
    ranks_pval <- c(GENE_A = 0.001, GENE_B = 0.01, GENE_A = 0.05)
    res_pval <- .select_best(ranks_pval, resolve_method = "rank")
    
    expect_equal(res_pval$keep_list$GENE_A, 1L) # Position 1 is highest ranked
    expect_equal(res_pval$keep_list$GENE_B, 2L) # B should be lower than A
    
    # If sorted by descending fold change, GENE_X at position 1 should be kept
    ranks_fc <- c(GENE_X = 5.0, GENE_Y = 3.0, GENE_X = 1.5)
    res_fc <- .select_best(ranks_fc, resolve_method = "rank")
    
    expect_equal(res_fc$keep_list$GENE_X, 1L) # Position 1 is highest ranked
})

test_that(".select_best() with resolve_method='rank' selects lowest rank", {
    ranks_test <- c(GENE_A = 25, GENE_B = 5, GENE_A = 3, GENE_B = 1, GENE_C = 10)
    ranks_test_srt <- sort(ranks_test)
    res_ranks <- .select_best(ranks_test_srt, resolve_method = "rank")
    
    expect_equal(res_ranks$keep_list$GENE_A, 2L) # A's lowest rank is pos 2, after sorting
    
    # Will produce a warning if ranks aren't sorted and "rank" is the selected method
    expect_warning(.select_best(ranks_test, resolve_method = "rank"),
                   "not monotonically sorted")
})

test_that(".select_best() with resolve_method='min' selects minimum signed value", {
    pval_test <- c(GENE_A = 0.05, GENE_B = 0.02, GENE_C = 0.01, GENE_C = 0.001)
    res <- .select_best(pval_test, resolve_method = "min")
    
    expect_equal(res$keep_list$GENE_C, 4L) # Only 2nd Gene C should be retained
    
    # Example with signed ranks
    ranks_multi <- c(GENE_X = 5, GENE_X = -5, GENE_Y = -3, GENE_X = 15)
    res_multi <- .select_best(ranks_multi, resolve_method = "min")
    
    expect_equal(res_multi$keep_list$GENE_X, 2L) # Index 2 has lowest value (5)
    expect_equal(res_multi$keep_list$GENE_Y, 3L) # Sanity check
})

test_that(".select_best() with resolve_method='max' selects maximum signed value", {
    pval_test <- c(GENE_A = 0.001, GENE_B = 0.02, GENE_A = 0.05)
    res <- .select_best(pval_test, resolve_method = "max")
    
    expect_equal(res$keep_list$GENE_A, 3L) # Index 3 has highest value (0.05)
    expect_equal(res$keep_list$GENE_B, 2L) # Only one occurrence
    
    # Example with signed ranks
    ranks_multi <- c(GENE_X = 5, GENE_X = -15, GENE_Y = 3, GENE_X = 10)
    res_multi <- .select_best(ranks_multi, resolve_method = "max")
    
    expect_equal(res_multi$keep_list$GENE_X, 4L) # Index 4 has max value (10)
    expect_equal(res_multi$keep_list$GENE_Y, 3L) # Sanity check
})

test_that(".select_best() resolve_method='rank' vs resolve_method='min' differ", {
    # Key distinction: 'rank' selects the entry closest to position 1 
    # (aka the first entry when sorted),
    # while 'min' selects the actual minimum value regardless of position.
    # This matters when a negative value appears at a later position than a 
    # positive value.
    
    # Sorted descending - typical for metrics like log fold change
    test_ranks <- c(GENE_X = 3, GENE_Z = 2, GENE_X = -1, GENE_Y = -2)
    
    res_rank <- .select_best(test_ranks, resolve_method = "rank")
    res_min  <- .select_best(test_ranks, resolve_method = "min")
    
    # 'rank' selects GENE_X at position 1
    expect_equal(res_rank$keep_list$GENE_X, 1L)
    expect_equal(test_ranks[res_rank$keep_list$GENE_X], c(GENE_X = 3))
    
    # 'min' selects GENE_X at position 3
    expect_equal(res_min$keep_list$GENE_X, 3L)
    expect_equal(test_ranks[res_min$keep_list$GENE_X], c(GENE_X = -1))
    
    # Confirm that results of 2 methods differ
    expect_false(res_rank$keep_list$GENE_X == res_min$keep_list$GENE_X)
    expect_false(identical(res_rank$keep_list, res_min$keep_list))
})

test_that(".select_best() handles pipe-delimited names", {
    ranks <- c("GENE_A|GENE_B" = 0.01, "GENE_C" = 0.02, "GENE_A" = 0.03)
    res <- .select_best(ranks)
    
    # .select_best() doesn't split names, they should remain in the output
    expect_true("GENE_A|GENE_B" %in% names(res$id_list))
    expect_equal(min(res$id_list$`GENE_A|GENE_B`), 1L)
    expect_equal(res$keep_list$`GENE_A|GENE_B`, 1L)
})

test_that(".select_best() preserves ranks vector structure when no pipe-delimited names", {
    simple_ranks <- c(A = 1, B = 2, C = 3)
    res <- .select_best(simple_ranks)
    
    expect_identical(res$ranks, simple_ranks)
})

test_that(".select_best() produces the expected number of results", {
    ranks <- abs(t_tests$log2_fc)
    names(ranks) <- t_tests$EntrezGeneSymbol
    ranks <- sort(ranks, decreasing = TRUE)
    res <- .select_best(ranks)
    
    n_names <- length(unique(names(ranks)))
    
    expect_length(res$keep_list, n_names)
})

test_that(".select_best() handles NA values in ranks", {
    na_ranks <- c(GENE_A = 0.01, GENE_A = NA, GENE_B = 0.05)
    res_min <- .select_best(na_ranks, resolve_method = "min")
    
    # which.min() with NA should still work, and is expected to
    # return the non-NA value
    expect_true("GENE_A" %in% names(res_min$keep_list))
    expect_true("GENE_B" %in% names(res_min$keep_list))
    expect_false(any(is.na(names(res_min$keep_list))))
})


####### .is.sorted() ########
test_that(".is.sorted() returns TRUE for ascending sorted vector", {
    expect_true(.is.sorted(c(1, 2, 3, 4, 5)))
})

test_that(".is.sorted() returns TRUE for descending sorted vector", {
    expect_true(.is.sorted(c(5, 4, 3, 2, 1)))
})

test_that(".is.sorted() returns FALSE for unsorted vector", {
    expect_false(.is.sorted(c(3, 1, 4, 1, 5)))
})

test_that(".is.sorted() returns TRUE or FALSE (never NA) when NAs are present", {
    na_asc  <- c(1, 2, NA, 4)   # ascending with NA
    na_desc <- c(4, NA, 2, 1)   # descending with NA
    na_unsorted <- c(3, NA, 1)  # unsorted with NA

    expect_true(!is.na(.is.sorted(na_asc)))
    expect_true(!is.na(.is.sorted(na_desc)))
    expect_true(!is.na(.is.sorted(na_unsorted)))
})
