# Setup ----

ranks <- t_tests$t_stat
names(ranks) <- t_tests$EntrezGeneSymbol
res <- .resolve_many_to_one(ranks = ranks, verbose = FALSE)

# Total unique groupings of genes & AptNames in the ex differential exp results.
# Does NOT include split heterodimers, that is tested by .splitIDs() and is
# not performed by .resolve_many_to_one()
genes_7k <- t_tests |>
    dplyr::select(AptName, EntrezGeneSymbol) |>
    dplyr::group_by(EntrezGeneSymbol) |>
    dplyr::summarize(n_apts = dplyr::n())

# Total number of genes targeted by >1 analyte
multimers <- dplyr::filter(genes_7k, n_apts > 1) |>
    dplyr::mutate(n_values_dropped = n_apts - 1)

# Reference values to test multimer filtering
n_multimapped_genes <- length(unique(multimers$EntrezGeneSymbol))
n_dropped_values <- sum(multimers$n_values_dropped)


# Testing -----
test_that("`.resolve_many_to_one()` produces the expected output with default args", {
    
    expect_type(res, "list")
    expect_named(res, c("results", "removed"))
    
    # Should be a named numeric vector, like the input
    expect_true(is.numeric(res$results))
    expect_named(res$results)
    
    # `removed` should be a named numeric vector of dropped values
    expect_true(is.numeric(res$removed))
    expect_named(res$removed)
    
    # Will differ for 5K/11K, the following tests apply to 7K only
    expect_length(res$removed, n_dropped_values)
    expect_length(res$results, nrow(genes_7k))
})

test_that("`.resolve_many_to_one()` selects the 'best' entry for each gene", {
    # Largest absolute value in ranks vector for KLK3 (the "best" should be kept)
    klk3_ranks <- ranks[which(names(ranks) == "KLK3")]
    klk3_top <- unname(klk3_ranks[which.max(abs(klk3_ranks))])
    klk3_kept <- res$results["KLK3"]
    
    # All non-top values should appear in removed
    klk3_removed <- res$removed[names(res$removed) == "KLK3"]
    
    # The kept value should match the entry with the largest absolute value
    expect_identical(unname(klk3_kept), klk3_top)
    expect_length(klk3_removed, length(klk3_ranks) - 1)
    expect_false(klk3_top %in% klk3_removed)
    
    # Same test as above, different gene that is seen in heterodimers
    # NOTE: Heterodimers aren't split by .resolve_many_to_one(), so actual CGA
    # heterodimers aren't tested here
    cga_ranks <- ranks[which(names(ranks) == "CGA")]
    cga_top <- unname(cga_ranks[which.max(abs(cga_ranks))])
    cga_kept <- res$results["CGA"]
    cga_removed <- res$removed[names(res$removed) == "CGA"]
    
    expect_identical(unname(cga_kept), cga_top)
    
    # All non-top values should appear in removed
    expect_length(cga_removed, length(cga_ranks) - 1)
    expect_false(cga_top %in% cga_removed)
})

test_that("`.resolve_many_to_one()` selects the expected entry, regardless of ranking metric", {
    ranks_fdr <- t_tests$fdr
    names(ranks_fdr) <- t_tests$EntrezGeneSymbol
    ranks_fdr <- sort(ranks_fdr, decreasing = TRUE)
    out_fdr <- .resolve_many_to_one(ranks = ranks_fdr, verbose = FALSE)
    
    klk3_vals <- ranks_fdr[which(names(ranks_fdr) == "KLK3")]
    top_klk3 <- max(klk3_vals)
    other_klk3 <- klk3_vals[klk3_vals != top_klk3]
    
    expect_identical(top_klk3, unname(out_fdr$results["KLK3"]))
    
    # All non-top values should be in removed
    klk3_removed <- out_fdr$removed[names(out_fdr$removed) == "KLK3"]
    expect_true(all(other_klk3 %in% klk3_removed))
})

# Snapshot test checking for correct numbers in the messages
test_that("`.resolve_many_to_one() produces the expected messaging`", {
    expect_snapshot(
        t <- .resolve_many_to_one(ranks = ranks, verbose = TRUE)
    )
})

test_that("`.resolve_many_to_one()` is silenced when `verbose = FALSE`", {
    expect_silent(
        .resolve_many_to_one(ranks = ranks, verbose = FALSE)
    )
})

test_that("`.resolve_many_to_one()` warns when ranks are unsorted and resolve_method = 'rank'", {
    idx <- withr::with_seed(505, sample(seq_len(length(ranks)), 
                                        size = length(ranks), 
                                        replace = FALSE))
    unsrt_ranks <- ranks[idx]
    
    expect_warning(
        .resolve_many_to_one(ranks = unsrt_ranks, resolve_method = "rank", verbose = FALSE), 
        "not monotonically sorted"
    )
})

test_that("`.resolve_many_to_one()` returns input unchanged when no duplicates present", {
    unique_ranks <- c(A = 10, B = 5, C = 3, D = 1)
    out <- .resolve_many_to_one(unique_ranks, verbose = FALSE)
    
    expect_identical(out$results, unique_ranks)
    expect_length(out$removed, 0)
})

test_that("`.resolve_many_to_one()` keeps instances of duplicate delimiter-separated genes", {
    # This function does not split heterodimers, so they should be retained!
    # Mock example ----------
    dup_ranks <- c(`A|A` = 10, A = 5, B = 5, C = 3, D = 1)
    out <- .resolve_many_to_one(dup_ranks, verbose = FALSE)
    
    expect_identical(out$results, dup_ranks)
    expect_length(out$removed, 0L) # No values are removed
    
    # Real-world example ----------
    df <- t_tests[grepl("GNAS", t_tests$EntrezGeneSymbol), ]
    dup_ranks_real <- df$t_stat
    names(dup_ranks_real) <- df$EntrezGeneSymbol
    out_real <- .resolve_many_to_one(dup_ranks_real, verbose = FALSE)

    expect_equal(out_real$results, dup_ranks_real)
    expect_length(out_real$removed, 0L)
    expect_identical(out_real$removed[1L], out_real$removed[2L])
})

test_that("`.resolve_many_to_one()` handles input when all entries are the same gene", {
    single_gene_ranks <- c(5, 3, 2, 1)
    names(single_gene_ranks) <- rep("GENE1", 4)
    
    out <- .resolve_many_to_one(single_gene_ranks, verbose = FALSE)
    
    # Should keep only the highest absolute value entry (5 in this case)
    expect_length(out$results, 1)
    expect_identical(unname(out$results), 5)
    expect_equal(names(out$results), "GENE1")
    
    # Should remove the other 3 entries
    expect_length(out$removed, 3)
    expect_true(all(out$removed %in% c(3, 2, 1)))
})

test_that("`.resolve_many_to_one()` produces different results based on 'resolve_method'", {
    # Case where "abs" and "rank" produce provably different results:
    # GENEA has three analytes: a large negative (-10) and two smaller positives (5, 3).
    # When sorted descending, "rank" keeps 5 (first occurrence).
    # "abs" always keeps -10 (largest absolute value), regardless of sort order.
    mock_ranks <- c(GENEA = 5, GENEA = 3, GENEA = -10, GENEB = 2, GENEC = 1)

    out_default <- .resolve_many_to_one(mock_ranks, verbose = FALSE)
    out_abs     <- .resolve_many_to_one(mock_ranks, resolve_method = "abs",
                                        verbose = FALSE)
    out_rank    <- suppressWarnings(
        .resolve_many_to_one(mock_ranks, 
                             resolve_method = "rank",
                             verbose = FALSE)
    )

    # Default should match "abs": keeps -10 for GENEA
    expect_identical(out_default$results, out_abs$results)

    # Default should NOT match "rank": "rank" keeps 5 for GENEA
    expect_false(identical(out_default$results, out_rank$results))

    # Explicit value check: "abs" keeps -10, "rank" keeps 5
    expect_equal(unname(out_default$results["GENEA"]), -10)
    expect_equal(unname(out_rank$results["GENEA"]), 5)
})

test_that("`.resolve_many_to_one()` 'results' and 'removed' are complementary, and reconstruct input", {
    # Every element in the input should appear in exactly one of results or removed
    mock_ranks <- c(GENEA = 5, GENEA = 3, GENEA = -10, GENEB = 2, GENEC = 1)
    out <- .resolve_many_to_one(mock_ranks, verbose = FALSE)

    reconstructed <- sort(c(out$results, out$removed))
    expect_identical(reconstructed, sort(mock_ranks))
})

test_that("`.resolve_many_to_one(resolve_method = 'min')` keeps minimum value per gene", {
    mock_ranks <- c(GENEA = 5, GENEA = -10, GENEA = 3, GENEB = 2)
    out <- .resolve_many_to_one(mock_ranks, resolve_method = "min", verbose = FALSE)

    expect_equal(unname(out$results["GENEA"]), -10)
    expect_true(all(c(5, 3) %in% out$removed))
})

test_that("`.resolve_many_to_one(resolve_method = 'max')` keeps maximum value per gene", {
    mock_ranks <- c(GENEA = 5, GENEA = -10, GENEA = 3, GENEB = 2)
    out <- .resolve_many_to_one(mock_ranks, resolve_method = "max", verbose = FALSE)

    expect_equal(unname(out$results["GENEA"]), 5)
    expect_true(all(c(-10, 3) %in% out$removed))
})

test_that("`.resolve_many_to_one()` handles NA values among duplicates", {
    na_ranks <- c(GENEA = NA_real_, GENEA = 5, GENEB = 2)
    out <- .resolve_many_to_one(na_ranks, resolve_method = "abs", verbose = FALSE)

    # "abs" should keep 5 (NA is not the largest absolute value)
    expect_equal(unname(out$results["GENEA"]), 5)
    expect_true(is.na(out$removed[["GENEA"]]))
})

test_that("`.resolve_many_to_one()` works regardless of input order", {
    asc_ranks <- sort(ranks, decreasing = FALSE)
    out <- .resolve_many_to_one(asc_ranks, verbose = FALSE)
    
    # Should work without error
    expect_type(out, "list")
    expect_named(out, c("results", "removed"))
    
    # With resolve_method = "abs" (default), sort order does not affect which
    # value is selected - the largest absolute value is always kept
    klk3_vals <- asc_ranks[names(asc_ranks) == "KLK3"]
    klk3_kept <- out$results["KLK3"]
    expect_identical(unname(klk3_kept), unname(klk3_vals[which.max(abs(klk3_vals))]))
})
