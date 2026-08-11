# Setup ----

meta       <- SomaDataIO::getAnalyteInfo(example_data_11k)
bg         <- meta$EntrezGeneSymbol
bg_apt     <- SomaDataIO::getAnalytes(example_data_11k)
deg        <- t_tests[t_tests$fc > 1 & t_tests$fdr < 0.05, ]$EntrezGeneSymbol
deg_apt    <- t_tests[t_tests$fc > 1 & t_tests$fdr < 0.05, ]$AptName

# Creating example pathways to use to test cust_paths arg
row_idx    <- withr::with_seed(101, sample(seq(1:nrow(pathway_map)), 75L, replace = FALSE))
cust_paths <- pathway_map[row_idx, ]
cust_paths$name <- "custom_path_1"
cust_paths <- df2list(cust_paths, name_col = "name", value_col = "gene_symbol")
cust_paths$custom_path_2 <- withr::with_seed(123, sample(
    pathway_map[pathway_map$group_code == "bp", ]$gene_symbol, 200L)
)

# Results object made w/ defaults
res_def    <- somaORA(features = deg, universe = bg)


# Testing -----
test_that("`somaORA()` produces the expected output with 'H' resource", {
    expect_snapshot(
        somaORA(features = deg, universe = bg, resource = "h", verbose = TRUE)
    )
})

# sanity check ---------
test_that("`somaORA()` has comparable results to `fgsea::fora()`", {
    # Save MSigDb Hallmark collection to use for example
    h_clxn <- pathway_map[pathway_map$group_code == "h", ]
    h_clxn <- df2list(h_clxn, name_col = "pathway_id", value_col = "gene_symbol")
    
    # Prep universe and background vectors
    u <- unique(t_tests$EntrezGeneSymbol)
    f <- u[1:50]
    
    # Perform ORA
    fora_res <- fgsea::fora(pathways = h_clxn,
                            genes = f,
                            universe = u)
    
    somaora_res <- somaORA(f, universe = u, 
                           cust_paths = h_clxn,
                           split_ids = FALSE)
    
    # Re-sort data to match result order
    fora_res <- as.data.frame(fora_res[order(fora_res$pathway), ])
    somaora_res <- somaora_res[order(somaora_res$pathway), ]
    
    # Filter to shared columns
    shared_cols <- intersect(colnames(fora_res), colnames(somaora_res))
    fora_res <- fora_res[, shared_cols]
    somaora_res <- somaora_res[, shared_cols]
    
    expect_identical(fora_res, somaora_res, 
                     ignore_attr = c("row.names", "waldo_opts")) # Need to ignore row names, since these were previously sorted
})

test_that("`somaORA()` will error if features aren't found in universe", {
    seqs <- c("seq.12345.67", "seq.891011.12", "seq.0000.00")
    expect_error(
        somaORA(features = c("NOTCH1", "HER2", "PIK3CA"), 
                universe = bg),
        "All elements in `features` must be found in `universe`."
    )
})

test_that("`somaORA(id_type = 'EntrezGeneID')` works with Entrez IDs", {
    # Build features and universe using Entrez gene IDs
    entrez_map <- unique(pathway_map[, c("gene_symbol", "entrez_id")])
    entrez_bg <- entrez_map$entrez_id[match(bg, entrez_map$gene_symbol)]
    entrez_bg <- unique(entrez_bg[!is.na(entrez_bg)])
    entrez_deg <- entrez_map$entrez_id[match(deg, entrez_map$gene_symbol)]
    entrez_deg <- entrez_deg[!is.na(entrez_deg)]
    
    res <- somaORA(features = entrez_deg, universe = entrez_bg,
                   id_type = "EntrezGeneID")
    
    expect_s3_class(res, "data.frame")
    expect_equal(colnames(res),
                 c("resource_code", "pathway_id", "pathway", "pval", "padj",
                   "foldEnrichment", "overlap", "size", "overlapFeatures"))
    expect_true(nrow(res) > 0L)
})

test_that("`somaORA(id_type = 'EntrezGeneID')` errors when gene symbols provided", {
    # Using gene symbol features with EntrezGeneID should find no overlap
    entrez_map <- unique(pathway_map[, c("gene_symbol", "entrez_id")])
    entrez_bg <- entrez_map$entrez_id[match(bg, entrez_map$gene_symbol)]
    entrez_bg <- unique(entrez_bg[!is.na(entrez_bg)])
    
    expect_error(
        somaORA(features = deg, universe = entrez_bg,
                id_type = "EntrezGeneID"),
        "All elements in `features` must be found in `universe`."
    )
})

# aptnames -------
test_that("`somaORA()` errors when only `features` are AptNames but `use_aptnames = FALSE`", {
    apts <- t_tests$AptName
    
    expect_error(
        somaORA(features = apts[1:20], universe = bg, resource = "h"),
        "All elements in `features` must be found in `universe`"
    )
})

test_that("`somaORA()` errors when AptNames are used but `use_aptnames = FALSE`", {
    apts <- t_tests$AptName

    expect_error(
        somaORA(features = apts[1:20], universe = bg_apt, resource = "h"),
        "`features` and `universe` appear to be SomaScan identifiers"
    )
})

test_that("`somaORA() errors when `use_aptnames = TRUE` but AptNames not provided" ,{
    apts <- SomaDataIO::getAnalytes(example_data_11k)
    
    expect_error(
        somaORA(features = deg, universe = bg, resource = "h", 
                use_aptnames = TRUE, col_meta_df = meta),
        "All features must be AptNames when `use_aptnames = TRUE`."
    )
})

# resource argument ----
test_that("`somaORA()` with resource = 'h' returns only Hallmark pathways", {
    h_pathways <- unique(path_name_lookup[path_name_lookup$pathway_id %in%
                             pathway_map[pathway_map$group_code == "h", ]$pathway_id, ]$pathway_name)
    res_h <- somaORA(features = deg, universe = bg, resource = "h")

    expect_true(all(res_h$resource_code == "h"))
    expect_true(all(res_h$pathway %in% h_pathways))
    expect_false(any(res_h$pathway %in%
                         unique(path_name_lookup[path_name_lookup$pathway_id %in%
                                    pathway_map[pathway_map$group_code == "bp", ]$pathway_id, ]$pathway_name)))
})

test_that("`somaORA()` has expected output with MSigDB resources", {
    res_c3  <- somaORA(features = deg, universe = bg, resource = "c3")
    c3_pths <- unique(path_name_lookup[path_name_lookup$pathway_id %in%
                         pathway_map[pathway_map$group_code == "c3", ]$pathway_id, ]$pathway_name)
    
    res_c8  <- somaORA(features = deg, universe = bg, resource = "c8")
    c8_pths <- unique(path_name_lookup[path_name_lookup$pathway_id %in%
                         pathway_map[pathway_map$group_code == "c8", ]$pathway_id, ]$pathway_name)
    
    expect_s3_class(res_c3, "data.frame")
    expect_true(all(res_c3$pathway %in% c3_pths))
    expect_true(all(res_c8$pathway %in% c8_pths))
})

test_that("`somaORA()` with `split_ids = TRUE` (default) splits heterodimers before gene set matching", {
    # Confirm heterodimers are present in the input features
    expect_true(any(grepl("\\|", deg)))

    res_split    <- somaORA(features = deg, universe = bg, resource = "h", split_ids = TRUE)
    res_no_split <- somaORA(features = deg, universe = bg, resource = "h", split_ids = FALSE)

    # Splitting heterodimers can only add more matchable features, not fewer
    expect_gte(sum(res_split$overlap,    na.rm = TRUE),
               sum(res_no_split$overlap, na.rm = TRUE))
})

test_that("`somaORA()` with `split_ids = FALSE` does not split heterodimers", {
    expect_true(any(grepl("\\|", deg)))

    # Standard gene sets don't contain "|" entries, so unsplit heterodimers
    # won't match any pathway gene and should not appear in overlapFeatures
    res <- somaORA(features = deg, universe = bg, resource = "h", split_ids = FALSE)
    expect_false(any(grepl("\\|", unlist(res$overlapFeatures))))
})

test_that("`somaORA()` runs without error when `unique_features = FALSE`", {
    # split_ids = TRUE can produce duplicate genes from split heterodimers;
    # unique_features = FALSE passes those duplicates through to fora()
    expect_true(any(grepl("\\|", deg)))
    expect_no_error(
        somaORA(features = deg, universe = bg, resource = "h",
                split_ids = TRUE, unique_features = FALSE)
    )
})

# cust_paths argument ----
test_that("`somaORA()` checks that `cust_paths` is a list", {
    path_vec <- unlist(cust_paths$custom_path_2)
    expect_error(
        somaORA(features = deg, universe = bg, cust_paths = path_vec),
        "`cust_paths` must be a list of named vectors."
    )
})

test_that("`somaORA()` errors when no elements in `features` are found in the selected gene sets", {
    apt_feats <- t_tests$AptName[1:20]
    apt_uni   <- t_tests$AptName

    expect_error(
        somaORA(features = apt_feats,
                universe = apt_uni,
                cust_paths = cust_paths),
        "No elements in `features` were found in the selected gene sets."
    )
})

test_that("`somaORA()` warns that `use_aptnames` is ignored when `cust_paths` is provided", {
    expect_warning(
        somaORA(features = deg,
                universe = bg,
                cust_paths = cust_paths,
                use_aptnames = TRUE),
        "`use_aptnames` is ignored when `cust_paths` is provided."
    )
})

test_that("`somaORA()` uses pathways provided to `cust_paths` argument", {
    res <- somaORA(features = deg, universe = bg, cust_paths = cust_paths)

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 2L)
    expect_true(all(res$pathway %in% c("custom_path_1", "custom_path_2")))
    # Custom pathways have no matching entry in pathway_map
    expect_true(all(is.na(res$pathway_id)))
    # overlapFeatures must be a subset of the corresponding custom pathway genes
    le1 <- res$overlapFeatures[[which(res$pathway == "custom_path_1")]]
    le2 <- res$overlapFeatures[[which(res$pathway == "custom_path_2")]]
    expect_true(all(le1 %in% cust_paths$custom_path_1))
    expect_true(all(le2 %in% cust_paths$custom_path_2))
})

test_that("`somaORA()` preferentially uses `cust_paths` over `resource`", {
    res <- somaORA(features = deg,
                   universe = bg,
                   resource = "h",
                   cust_paths = cust_paths)

    # Should only contain custom pathways, not Hallmark pathways
    expect_equal(nrow(res), 2L)
    expect_true(all(res$pathway %in% c("custom_path_1", "custom_path_2")))
    expect_false(any(grepl("^HALLMARK", res$pathway)))
})

test_that("`somaORA()` works with `use_aptnames = TRUE` and AptName features", {
    apt_features <- bg_apt[1:100]

    res_apt <- somaORA(features = apt_features,
                       universe = bg_apt,
                       resource = "h",
                       use_aptnames = TRUE,
                       col_meta_df = meta)

    expect_s3_class(res_apt, "data.frame")
    expect_equal(colnames(res_apt),
                 c("resource_code", "pathway_id", "pathway", "pval", "padj",
                   "foldEnrichment", "overlap", "size", "overlapFeatures"))
    expect_true(nrow(res_apt) > 0L)
    # overlapFeatures should all be AptNames
    expect_true(all(SomaDataIO::is.AptName(unlist(res_apt$overlapFeatures))))
})

test_that("`somaORA()` errors when `use_aptnames = TRUE` but `col_meta_df` not provided", {
    apt_features <- bg_apt[1:100]
    expect_error(
        somaORA(features = apt_features, universe = bg_apt,
                resource = "h", use_aptnames = TRUE),
        "`col_meta_df` must be provided when `use_aptnames = TRUE`."
    )
})

test_that("`somaORA()` returns fewer results with a higher `min_feats` threshold", {
    res_low  <- somaORA(features = deg, universe = bg, resource = "h", min_feats = 5L)
    res_high <- somaORA(features = deg, universe = bg, resource = "h", min_feats = 50L)

    expect_gt(nrow(res_low), nrow(res_high))
    # Every retained pathway must have at least min_feats members post-filtering
    expect_true(all(res_high$size >= 50L))
})

test_that("`somaORA()` returns fewer results with a lower `max_feats` threshold", {
    res_high <- somaORA(features = deg, universe = bg, resource = "h", max_feats = 500L)
    res_low  <- somaORA(features = deg, universe = bg, resource = "h", max_feats = 50L)

    expect_gt(nrow(res_high), nrow(res_low))
    # Every retained pathway must have no more than max_feats members post-filtering
    expect_true(all(res_low$size <= 50L))
})

# verbosity and messaging -------
test_that("`somaORA()` is appropriately silenced when `verbose = FALSE`", {
    expect_no_message(
        somaORA(features = deg, universe = bg, resource = "h", verbose = FALSE)
    )
})

test_that("`somaORA()` produces a progress message when `verbose = TRUE`", {
    expect_snapshot(
       t <- somaORA(features = deg, universe = bg, resource = "h", verbose = TRUE)
    )
})

test_that("`somaORA()` results are sorted by p-adj (ascending)", {
    expect_true(all(diff(res_def$padj) >= 0))
})

test_that("`somaORA()` handles duplicate elements in universe", {
    dup_bg <- c(bg, bg[1:10])
    expect_no_error(
        somaORA(features = deg, universe = dup_bg, resource = "h")
    )
    # Output should be the same, since duplicate universe elements are always removed
    res_dup    <- somaORA(features = deg,
                          universe = dup_bg, 
                          resource = "h")
    res_no_dup <- somaORA(features = deg,
                          universe = bg,
                          resource = "h")
    expect_identical(res_dup$pval, res_no_dup$pval)
})

test_that("`somaORA()` removes NA values from `universe` without error", {
    # NAs in universe should be silently dropped (same behavior as empty strings)
    bg_with_na <- c(bg, NA_character_, NA_character_)
    
    expect_no_error(
        somaORA(features = deg, universe = bg_with_na, resource = "h")
    )
    res_na  <- somaORA(features = deg, universe = bg_with_na, resource = "h")
    res_clean <- somaORA(features = deg, universe = bg, resource = "h")
    
    expect_identical(res_na$pval, res_clean$pval)
})

test_that("`somaORA()` errors gracefully when `features` is empty", {
    expect_error(
        somaORA(features = character(0), universe = bg, resource = "h"),
        "`features` must contain at least one element."
    )
})

