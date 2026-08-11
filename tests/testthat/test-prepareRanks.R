# Setup ---------

# "Real-world" example
sorted_data <- t_tests[order(t_tests$t_stat), ]
sorted_vec <- sorted_data$t_stat
names(sorted_vec) <- sorted_data$EntrezGeneSymbol

noDimer_data <- sorted_data[!grepl("\\|", sorted_data$EntrezGeneSymbol), ]
noDimer_vec <- noDimer_data$t_stat
names(noDimer_vec) <- noDimer_data$EntrezGeneSymbol

# Simple mock input with 1:1 mapping and no heterodimers
sim_stats    <- c(1.0, 2.0, 3.0)
sim_features <- c("GeneA", "GeneB", "GeneC")

# Mock input with pipe-delimited heterodimers
het_stats    <- c(1.5, 2.0, 3.0)
het_features <- c("GeneA|GeneB", "GeneC", "GeneD")

# Mock input with multimapping (but no heterodimers)
dup_stats    <- c(1.0, 2.0, 3.0, 4.0)
dup_features <- c("GeneA", "GeneB", "GeneA", "GeneC")


# Testing --------
# This is just a wrapper function, goal of tests is to make sure this
# function doesn't add anything unexpected. All internal functions have their
# own testing suites

test_that("`prepareRanks()` returns expected output with example data", {
    ranks <- t_tests$t_stat
    names(ranks) <- t_tests$EntrezGeneSymbol
    ranks <- sort(ranks)
    
    start_len <- length(ranks)
    split_ranks <- .splitIDs(ranks)
    split_len <- length(split_ranks)
    filt_ranks <- .resolve_many_to_one(split_ranks, verbose = FALSE)
    filt_len <- length(filt_ranks$results)
    rmv_len <- length(filt_ranks$removed)
    
    prepped_ranks <- prepareRanks(stats = sorted_data$t_stat,
                                  features = sorted_data$EntrezGeneSymbol)
    
    # Should be same length as internal steps run individually
    expect_length(prepped_ranks, split_len - rmv_len)
})

test_that("`prepareRanks()` returns a named numeric vector", {
    result <- prepareRanks(stats = sorted_data$t_stat,
                           features = sorted_data$EntrezGeneSymbol,
                           verbose = FALSE)
    expect_true(is.numeric(result))
    expect_named(result)
})

test_that("`prepareRanks()` returns input vector when no heterodimers or duplication present", {
    result <- prepareRanks(stats = sim_stats,
                           features = sim_features,
                           verbose = FALSE)
    named_vec <- sim_stats
    names(named_vec) <- sim_features
    
    expect_equal(named_vec, result)
})

test_that("`prepareRanks()` output never has '|' in names when 'split_ids' = TRUE", {
    result <- prepareRanks(stats = sorted_data$t_stat,
                           features = sorted_data$EntrezGeneSymbol,
                           verbose = FALSE)
    expect_false(any(grepl("\\|", names(result))))
})

test_that("`prepareRanks()` output has no duplicate names", {
    result <- prepareRanks(stats = sorted_data$t_stat,
                           features = sorted_data$EntrezGeneSymbol,
                           verbose = FALSE)
    expect_false(any(duplicated(names(result))))
})

test_that("`prepareRanks(split_ids = TRUE)` expands heterodimers", {
    result <- prepareRanks(stats = het_stats,
                           features = het_features,
                           split_ids = TRUE,
                           resolve_multimapping = FALSE,
                           verbose = FALSE)
    # "GeneA|GeneB" -> "GeneA" and "GeneB", so 4 elements total
    expect_length(result, 4L)
    expect_true(all(c("GeneA", "GeneB", "GeneC", "GeneD") %in% names(result)))
})

test_that("`prepareRanks(split_ids = FALSE)` preserves heterodimer names", {
    result <- prepareRanks(stats = het_stats,
                           features = het_features,
                           split_ids = FALSE,
                           resolve_multimapping = FALSE,
                           verbose = FALSE)
    expect_length(result, 3L)
    expect_true("GeneA|GeneB" %in% names(result))
})


test_that("`prepareRanks(resolve_multimers = FALSE)` errors when duplicates exist", {
    expect_error(
        prepareRanks(stats = dup_stats,
                     features = dup_features,
                     split_ids = FALSE,
                     resolve_multimapping = FALSE,
                     verbose = FALSE),
        "Duplicate values found"
    )
})


test_that("`prepareRanks()` returns early when no duplicates remain after splitting", {
    result <- prepareRanks(stats = het_stats,
                           features = het_features,
                           split_ids = TRUE,
                           resolve_multimapping = TRUE,
                           verbose = FALSE)
    # After splitting, all names are unique -> early return, no filtering needed
    expect_length(result, 4L)
    expect_false(any(duplicated(names(result))))
})

test_that("`prepareRanks()` returns early with simple unique input", {
    result <- prepareRanks(stats = sim_stats,
                           features = sim_features,
                           verbose = FALSE)
    expect_identical(unname(result), sim_stats)
    expect_identical(names(result), sim_features)
})

test_that("`prepareRanks(verbose = TRUE)` emits expected messages", {
    expect_snapshot(
        t <- prepareRanks(stats = dup_stats,
                          features = dup_features,
                          split_ids = FALSE,
                          resolve_multimapping = TRUE,
                          verbose = TRUE)
    )
})

test_that("`prepareRanks(verbose = FALSE)` is silent", {
    expect_silent(
        prepareRanks(stats = sorted_data$t_stat,
                     features = sorted_data$EntrezGeneSymbol,
                     verbose = FALSE)
    )
})

test_that("`prepareRanks()` output identical to ranking vector in somaPrGSEA() output", {
    pr_res <- prepareRanks(stats = noDimer_data$t_stat,
                           features = noDimer_data$EntrezGeneSymbol,
                           verbose = FALSE)
    gp_res <- somaPrGSEA(ranks = noDimer_vec,
                         verbose = FALSE)
    
    expect_equal(pr_res, gp_res$final_ranks)
})

test_that("`prepareRanks()` warns about NAs, with the correct count", {
    na_stats    <- c(NA, 2.0, NA, 4.0)
    na_features <- c("GeneA", "GeneB", "GeneC", "GeneD")

    # Snapshots are easier than trying to capture the formatting used by cli
    expect_snapshot(
        prepareRanks(stats = na_stats, features = na_features, verbose = FALSE),
    )
    
    # Should only report 1 NA, with no pluralization
    expect_snapshot(
        prepareRanks(stats = na_stats[-1], features = na_features[-1], verbose = FALSE),
    )
})

test_that("`prepareRanks()` errors when `stats` and `features` have different lengths", {
    expect_error(
        prepareRanks(stats = c(1.0, 2.0, 3.0),
                     features = c("GeneA", "GeneB"),
                     verbose = FALSE),
        "`stats` and `features` must have the same length."
    )
    expect_error(
        prepareRanks(stats = c(1.0),
                     features = character(0),
                     verbose = FALSE),
        "`stats` and `features` must have the same length."
    )
})

test_that("`prepareRanks()` works correctly with all `resolve_method` options", {
    # `resolve_method = 'rank'` requires sorted input.
    # For the others, any order is fine.
    
    # --- 'rank': keeps first occurrence in a sorted vector ---
    # Sort descending: 5.0, 3.0, 1.0 -> GeneA at pos 1 (5.0) is kept
    stats_sorted    <- c(5.0, 3.0, 1.0)
    features_sorted <- c("GeneA", "GeneA", "GeneB")
    res_rank <- prepareRanks(stats_sorted, features_sorted,
                             split_ids = FALSE, resolve_multimapping = TRUE,
                             resolve_method = "rank", verbose = FALSE)
    expect_equal(unname(res_rank["GeneA"]), 5.0)  # first (highest rank) kept
    
    # --- 'min': keeps the minimum numeric value for each gene ---
    stats_multi    <- c(5.0, 1.0, 3.0)
    features_multi <- c("GeneA", "GeneB", "GeneA")
    res_min <- prepareRanks(stats_multi, features_multi,
                            split_ids = FALSE, resolve_multimapping = TRUE,
                            resolve_method = "min", verbose = FALSE)
    # GeneA has values 5.0 and 3.0; min within the group selects 3.0
    expect_equal(unname(res_min["GeneA"]), 3.0)
    expect_equal(unname(res_min["GeneB"]), 1.0)  # GeneB has only one value
    
    # --- 'max': keeps the maximum numeric value ---
    res_max <- prepareRanks(stats_multi, features_multi,
                            split_ids = FALSE, resolve_multimapping = TRUE,
                            resolve_method = "max", verbose = FALSE)
    expect_equal(unname(res_max["GeneA"]), 5.0)
    
    # --- 'abs': keeps the entry with the largest absolute value ---
    # GeneA has -5.0 and 1.0 -> abs selects -5.0
    stats_abs    <- c(-5.0, 1.0, 1.0)
    features_abs <- c("GeneA", "GeneB", "GeneA")
    res_abs <- prepareRanks(stats_abs, features_abs,
                            split_ids = FALSE, resolve_multimapping = TRUE,
                            resolve_method = "abs", verbose = FALSE)
    expect_equal(unname(res_abs["GeneA"]), -5.0)
})

# Output sorting tests --------
test_that("`prepareRanks(resolve_method = 'rank')` output preserves input sort direction", {
    # Ascending input -> ascending output
    stats_asc    <- c(1.0, 2.0, 3.0, 5.0, 7.0)
    features_asc <- c("GeneA", "GeneB", "GeneA", "GeneC", "GeneD")
    res_asc <- prepareRanks(stats_asc, features_asc,
                            split_ids = FALSE, resolve_multimapping = TRUE,
                            resolve_method = "rank", verbose = FALSE)
    expect_equal(sort(res_asc), res_asc)

    # Descending input -> descending output
    stats_desc    <- c(7.0, 5.0, 3.0, 2.0, 1.0)
    features_desc <- c("GeneA", "GeneB", "GeneA", "GeneC", "GeneD")
    res_desc <- prepareRanks(stats_desc, features_desc,
                             split_ids = FALSE, resolve_multimapping = TRUE,
                             resolve_method = "rank", verbose = FALSE)
    expect_equal(sort(res_desc, decreasing = TRUE), res_desc)
})

test_that("`prepareRanks()` output passes to `somaPrGSEA()` without warnings or errors", {
    result <- prepareRanks(stats = sorted_data$t_stat,
                           features = sorted_data$EntrezGeneSymbol,
                           split_ids = FALSE,
                           verbose = FALSE)
    expect_no_condition(
        somaPrGSEA(ranks = result,
                   resource = "h",
                   split_ids = FALSE,
                   resolve_multimapping = FALSE,
                   verbose = FALSE)
    )
})

