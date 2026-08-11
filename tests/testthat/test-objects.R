
test_that("`pathway_map` object has the expected characteristics", {
    expect_s3_class(pathway_map, "data.frame")
    expect_named(pathway_map, c("gene_symbol", "entrez_id", "pathway_id", "group_code"))
    expect_s3_class(pathway_map$group_code, "factor")
    
    # Row counts may shift slightly across upstream database releases.
    # Baselines were captured from msigdbr 26.1.1 / GO.db 3.23.0.
    # Tolerances of ±20% flag gross structural problems without requiring
    # exact counts to be maintained after each re-generation.
    baseline_total  <- 4489940L
    baseline_bp     <- 1284385L
    baseline_mf     <-  270956L
    baseline_other  <- 3072275L
    tol <- 0.20

    expect_equal(ncol(pathway_map), 4L)
    expect_gte(nrow(pathway_map), baseline_total * (1 - tol))
    expect_lte(nrow(pathway_map), baseline_total * (1 + tol))

    n_bp    <- nrow(pathway_map[pathway_map$group_code == "bp", ])
    n_mf    <- nrow(pathway_map[pathway_map$group_code == "mf", ])
    n_other <- nrow(pathway_map[!pathway_map$group_code %in% c("bp", "mf"), ])
    expect_gte(n_bp,    baseline_bp    * (1 - tol))
    expect_lte(n_bp,    baseline_bp    * (1 + tol))
    expect_gte(n_mf,    baseline_mf    * (1 - tol))
    expect_lte(n_mf,    baseline_mf    * (1 + tol))
    expect_gte(n_other, baseline_other * (1 - tol))
    expect_lte(n_other, baseline_other * (1 + tol))
})

test_that("`pathway_map` object contains expected GO/MSigDB subcollections", {
    colxns <- table(pathway_map$group_code)
    expect_named(colxns, c("bp", "c1", "c2", "c3", "c4", "c6", "c7", "c8", "h", "mf"))
})

test_that("`group_code_lookup` internal object has expected structure", {
    expect_s3_class(group_code_lookup, "data.frame")
    expect_equal(nrow(group_code_lookup), 10L)
    expect_named(group_code_lookup, c("resource", "collection", "group_code", "label"))
    expect_setequal(group_code_lookup$group_code,
                    c("bp", "mf", "h", "c1", "c2", "c3", "c4", "c6", "c7", "c8"))
})

test_that("`path_name_lookup` internal object has expected structure", {
    expect_s3_class(path_name_lookup, "data.frame")
    expect_named(path_name_lookup, c("pathway_id", "pathway_name"))
    expect_equal(nrow(path_name_lookup), length(unique(pathway_map$pathway_id)))
})
