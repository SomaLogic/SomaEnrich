# Setup ----
meta            <- getAnalyteInfo(example_data_11k)
ranks           <- t_tests$t_stat
names(ranks)    <- t_tests$EntrezGeneSymbol
ranks           <- sort(ranks, decreasing = TRUE)
ranks_noDimer    <- ranks[!grepl("\\|", names(ranks))]
row_idx         <- withr::with_seed(101, sample(seq(1:nrow(pathway_map)), 75L,
                                                replace = FALSE))
cust_paths      <- pathway_map[row_idx, ]
cust_paths$name <- "custom_path_1"
cust_paths      <- df2list(cust_paths, name_col = "name", value_col = "gene_symbol")
cust_paths$custom_path_2 <- withr::with_seed(101, sample(
    pathway_map[pathway_map$group_code == "bp", ]$gene_symbol, 200L)
)

# Will produce ties warning - suppress for this step
res_def <- suppressWarnings(somaPrGSEA(ranks = ranks))


# Testing -----
test_that("`somaPrGSEA()` output has expected structure", {
  expect_type(res_def, "list")
  expect_named(res_def, c("results", "final_ranks"))
  expect_s3_class(res_def$results, "data.frame")
  expect_equal(colnames(res_def$results),
               c("resource_code", "pathway_id", "pathway", "pval", "padj", 
                 "log2err", "ES", "NES", "leadingEdge", 
                 "starting_set_size", "final_set_size"))
  # Final_ranks is a named numeric vector with no duplicate names
  expect_false(is.null(names(res_def$final_ranks)))
  expect_false(any(duplicated(names(res_def$final_ranks))))
})

test_that("`somaPrGSEA()` has expected results with 'H' resource", {
  # Use noDimer vector to avoid any fgsea-created warning messages
  res <- somaPrGSEA(ranks = ranks_noDimer, 
                    resource = "h", 
                    verbose = FALSE,
                    nproc = 1) # Necessary for perfect reproducibility
    
  expect_snapshot(res$results)
  expect_snapshot(res$final_ranks[1:1000])
})

test_that("`somaPrGSEA(resource = 'h')` returns only Hallmark pathways", {
    h_pathways <- unique(path_name_lookup[path_name_lookup$pathway_id %in%
                             pathway_map[pathway_map$group_code == "h", ]$pathway_id, ]$pathway_name)
    res_h <- somaPrGSEA(ranks = ranks_noDimer, resource = "h")
    
    expect_true(all(res_h$results$resource_code == "h"))
    expect_true(all(res_h$results$pathway %in% h_pathways))
    expect_false(any(res_h$results$pathway %in%
                         unique(path_name_lookup[path_name_lookup$pathway_id %in%
                                    pathway_map[pathway_map$group_code == "bp", ]$pathway_id, ]$pathway_name)))
})

# sanity check ---------
test_that("`somaPrGSEA()` has comparable results to fgsea", {
    
    # Save MSigDb Hallmark collection to use for example
    h_clxn <- pathway_map[pathway_map$group_code == "h", ]
    h_clxn <- df2list(h_clxn, name_col = "pathway_id", value_col = "gene_symbol")
    
    # Prep ranked stats
    r <- t_tests$t_stat
    names(r) <- t_tests$EntrezGeneSymbol
    r <- r[!duplicated(names(r))] # Arbitrarily taking the 1st occurrence of all duplicates
    r <- sort(r, decreasing = TRUE)
    
    # Perform GSEA
    fgsea_res <- withr::with_seed(42, fgsea::fgsea(pathways = h_clxn, stats = r))
    somagsea_res <- somaPrGSEA(r, cust_paths = h_clxn,
                               resolve_multimapping = FALSE,
                               split_ids = FALSE)
    
    # Re-sort data to match result order
    fgsea_res <- as.data.frame(fgsea_res[order(fgsea_res$pathway), ])
    somagsea_res$results <- somagsea_res$results[order(somagsea_res$results$pathway), ]
    
    # Filter to shared columns
    shared_cols <- intersect(colnames(fgsea_res), colnames(somagsea_res$results))
    fgsea_res <- fgsea_res[, shared_cols]
    somagsea_res$results <- somagsea_res$results[, shared_cols]
    
    expect_identical(somagsea_res$final_ranks, r)
    expect_identical(fgsea_res, somagsea_res$results, 
                     ignore_attr = c("row.names", "waldo_opts")) # Need to ignore row names, since these were previously sorted
})

# error catching ------
test_that("`somaPrGSEA()` errors when IDs don't match in pathways and ranks", {
    err_ranks <- t_tests$t_stat
    names(err_ranks) <- t_tests$UniProt
    err_ranks <- sort(err_ranks, decreasing = TRUE)
    
    expect_error(
        somaPrGSEA(ranks = err_ranks, resource = "h"),
        "No elements in `ranks` found in the selected gene sets."
    )
})

test_that("`somaPrGSEA()` warns when the ranks are unsorted AND multimers are present", {
    unsrt_ranks <- withr::with_seed(123, sample(ranks_noDimer[1:15], size = 10,
                                                replace = TRUE))
    # Sanity check
    # .is.sorted(unsrt_ranks)
    expect_warning(
        somaPrGSEA(ranks = unsrt_ranks, resolve_method = "rank"), 
        "not monotonically sorted"
    )
})

test_that("`somaPrGSEA()` errors when no elements in `ranks` are found in the chosen resource", {
    ranks <- seq(1:5L)
    names(ranks) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
    
    expect_error(
        somaPrGSEA(ranks = ranks, resource = "h"),
        "No elements in `ranks` found in the selected gene sets."
    )
})

test_that("`somaPrGSEA()` checks that `cust_paths` is a list", {
    path_vec <- unlist(cust_paths$custom_path_2)
    expect_error(
        somaPrGSEA(ranks = ranks_noDimer,
                   cust_paths = path_vec),
        "`cust_paths` must be a list of named vectors."
    )
})

test_that("`somaPrGSEA()` warns that `use_aptnames` is ignored when `cust_paths` is provided", {
    expect_warning(
        somaPrGSEA(ranks = ranks_noDimer,
                   cust_paths = cust_paths,
                   use_aptnames = TRUE,
                   split_ids = FALSE),
        "`use_aptnames` is ignored when `cust_paths` is provided."
    )
})

# verbosity -------
test_that("`somaPrGSEA(verbose = TRUE)` emits expected messages", {
    # Use input with no heterodimers to avoid testing fgsea warnings
    expect_snapshot(
        t <- somaPrGSEA(ranks = ranks_noDimer,
                        verbose = TRUE)
    )
})

test_that("`somaPrGSEA()` is appropriately silenced when `verbose = FALSE`", {
    expect_no_error(
        somaPrGSEA(ranks = ranks_noDimer, 
                   resource = "h",
                   verbose = FALSE,
                   split_ids = FALSE)
    )
})

# custom gene sets/pathways -----
test_that("`somaPrGSEA()` uses dataset provided to `cust_paths` argument", {
    res <- somaPrGSEA(ranks = ranks_noDimer,
                      cust_paths = cust_paths,
                      split_ids = FALSE)
    le1 <- res$leadingEdge[[1L]]
    le2 <- res$leadingEdge[[2L]]
    
    # Should be same length as cust_paths list
    expect_equal(nrow(res$results), 2L)
    
    # Shouldn't be the "bp" default, should be "custom" x number of gene sets
    expect_equal(res$results$resource_code, c("custom", "custom"))
    expect_true(all(le1 %in% cust_paths$custom_path_1))
    expect_true(all(le2 %in% cust_paths$custom_path_2))
})

test_that("`somaPrGSEA()` preferentially uses `cust_paths` over `resource`", {
    res <- somaPrGSEA(ranks = ranks_noDimer, 
                      resource = "h", 
                      cust_paths = cust_paths,
                      split_ids = FALSE)
    
    # Should NOT contain pathways from GO Molecular Function
    expect_equal(nrow(res$results), 2L)
    expect_equal(res$results$resource_code, c("custom", "custom"))
    expect_true(all(res$results$pathway %in% c("custom_path_1", "custom_path_2")))
})

test_that("`somaPrGSEA()` warns about duplicates when `resolve_multimapping = FALSE`", {
    expect_error(
        somaPrGSEA(ranks = ranks_noDimer, resolve_multimapping = FALSE),
        "Please remove duplicates before proceeding"
        )
})

test_that("`somaPrGSEA(id_type = 'EntrezGeneID')` works with Entrez IDs", {
    # Build ranks named with Entrez gene IDs
    t_split <- tidyr::separate_rows(t_tests, EntrezGeneSymbol, sep = "\\|")
    entrez_lookup <- pathway_map$entrez_id[match(t_split$EntrezGeneSymbol,
                                                 pathway_map$gene_symbol)]
    entrez_ranks <- t_split$t_stat
    names(entrez_ranks) <- entrez_lookup
    entrez_ranks <- entrez_ranks[!is.na(names(entrez_ranks))]
    entrez_ranks <- sort(entrez_ranks, decreasing = TRUE)
    
    # Will produce warning about ties
    res <- suppressWarnings(somaPrGSEA(ranks = entrez_ranks, 
                                       id_type = "EntrezGeneID",
                                       resource = "h"))
    
    expect_s3_class(res$results, "data.frame")
    expect_equal(colnames(res$results),
                 c("resource_code", "pathway_id", "pathway", "pval", "padj",
                   "log2err", "ES", "NES", "leadingEdge",
                   "starting_set_size", "final_set_size"))
    expect_true(nrow(res$results) > 0L)
    expect_true(all(grepl("^HALLMARK", res$results$pathway)))
})

test_that("`somaPrGSEA(id_type = 'EntrezGeneID')` errors when gene symbols provided", {
    # Using gene symbol-named ranks with EntrezGeneID should fail
    expect_error(
        somaPrGSEA(ranks = ranks_noDimer, id_type = "EntrezGeneID", resource = "h"),
        "No elements in `ranks` found in the selected gene sets."
    )
})

# use_aptnames argument ----
test_that("`somaPrGSEA()` works with use_aptnames = TRUE and AptName-named ranks", {
    ranks_apt <- t_tests$t_stat
    names(ranks_apt) <- t_tests$AptName
    ranks_apt <- sort(ranks_apt, decreasing = TRUE)

    res_apt <- somaPrGSEA(ranks = ranks_apt,
                          col_meta_df = meta,
                          use_aptnames = TRUE,
                          resource = "h")

    expect_s3_class(res_apt$results, "data.frame")
    expect_true(all(res_apt$results$resource_code == "h"))
    # final_ranks should be in AptName (seq.NNNN.NN) format
    expect_true(all(SomaDataIO::is.AptName(names(res_apt$final_ranks))))
    expect_true((all(SomaDataIO::is.AptName(res_apt$results$leadingEdge[[1]]))))
})

test_that("`somaPrGSEA()` errors when AptNames are in ranks but `use_aptnames = FALSE`", {
    ranks_apt <- t_tests$t_stat
    names(ranks_apt) <- t_tests$AptName
    ranks_apt <- sort(ranks_apt, decreasing = TRUE)

    expect_error(
        somaPrGSEA(ranks = ranks_apt, resource = "h"),
        "`ranks` appears to be SomaScan identifiers, did you set `use_aptnames = TRUE`?"
    )
})

test_that("`somaPrGSEA()` errors when use_aptnames = TRUE but `col_meta_df` not provided", {
    ranks_apt <- t_tests$t_stat
    names(ranks_apt) <- t_tests$AptName
    ranks_apt <- sort(ranks_apt, decreasing = TRUE)

    expect_error(
        somaPrGSEA(ranks = ranks_apt, use_aptnames = TRUE, resource = "h"),
        "`col_meta_df` must be provided when `use_aptnames = TRUE`."
    )
})

# min_feats / max_feats arguments ----
test_that("`somaPrGSEA()` returns fewer results with a higher `min_feats` threshold", {
    res_low  <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", min_feats = 5L)
    res_high <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", min_feats = 50L)

    expect_gt(nrow(res_low$results), nrow(res_high$results))
    # Every retained pathway must have at least min_feats members post-filtering
    expect_true(all(res_high$results$final_set_size >= 50L))
})

test_that("`somaPrGSEA()` returns fewer results with a lower `max_feats` threshold", {
    res_high <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", max_feats = 500L)
    res_low  <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", max_feats = 50L)

    expect_gt(nrow(res_high$results), nrow(res_low$results))
    # Every retained pathway must have no more than max_feats members post-filtering
    expect_true(all(res_low$results$final_set_size <= 50L))
})

# split_ids argument ----
test_that("`somaPrGSEA()` with `split_ids = FALSE` retains heterodimers in `final_ranks`", {
    # Confirm heterodimers are present in the input ranks
    expect_true(any(grepl("\\|", names(ranks))))

    res <- somaPrGSEA(ranks = ranks, resource = "h", split_ids = FALSE)
    expect_true(any(grepl("\\|", names(res$final_ranks))))
})

test_that("`somaPrGSEA(split_ids = TRUE)` removes heterodimers from `final_ranks`", {
    # Confirm heterodimers are present in the input ranks before processing
    expect_true(any(grepl("\\|", names(ranks))))
    expect_warning(somaPrGSEA(ranks = ranks, resource = "h"),
                   "There are ties in the preranked stats")
    expect_false(any(grepl("\\|", names(res_def$final_ranks))))
})

# resolve_multimapping argument ----
test_that("`somaPrGSEA(resolve_multimapping = FALSE)` passes pre-unique ranks through unchanged", {
    # Remove any duplicated gene names so no processing is needed
    unique_ranks <- ranks_noDimer[!duplicated(names(ranks_noDimer))]

    res <- somaPrGSEA(ranks = unique_ranks,
                      resource = "h",
                      split_ids = FALSE,
                      resolve_multimapping = FALSE)

    expect_equal(res$final_ranks, unique_ranks)
})

test_that("`somaPrGSEA(resolve_method = 'min')` produces different final_ranks than default 'abs'", {
    # 'rank' takes the first (highest) value per gene in the sorted vector;
    # 'min' takes the actual minimum value, which may be negative or smaller
    # than the first occurrence. For a descending t-stat vector with
    # multi-mapped genes having distinct values, these will differ.
    res_rank <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", resolve_method = "rank")
    res_min  <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", resolve_method = "min")

    expect_false(identical(res_rank$final_ranks, res_min$final_ranks))
})

test_that("`somaPrGSEA()` produces identical results when the same seed is used", {
    res1 <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", seed = 42L)
    res2 <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", seed = 42L)

    expect_identical(res1$results$pval, res2$results$pval)
    expect_identical(res1$results$NES,  res2$results$NES)
})

test_that("`somaPrGSEA()` produces different results when a different seed is used", {
    res_s42 <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", seed = 42L)
    res_s99 <- somaPrGSEA(ranks = ranks_noDimer, resource = "h", seed = 99L)

    # Different seeds should produce different p-values for the stochastic algorithm
    expect_false(identical(res_s42$results$pval, res_s99$results$pval))
})

test_that("`somaPrGSEA(score_type = 'pos')` handles all-positive ranks without error", {
    # Create a unique, all-positive ranks vector (no duplicates)
    pos_ranks <- abs(ranks_noDimer)
    pos_ranks <- pos_ranks[!duplicated(names(pos_ranks))]
    pos_ranks <- sort(pos_ranks, decreasing = TRUE)

    expect_no_error(
        somaPrGSEA(ranks = pos_ranks, resource = "h",
                   split_ids = FALSE, resolve_multimapping = FALSE)
    )
})

