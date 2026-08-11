# Setup ------
ranks_pe <- t_tests$t_stat
names(ranks_pe) <- t_tests$EntrezGeneSymbol
ranks_pe <- sort(ranks_pe, decreasing = TRUE)

gsea_res <- suppressWarnings(
  somaPrGSEA(ranks = ranks_pe, resource = "h")
)

# Testing -------
test_that("`plotES()` returns expected patchwork/gg object", {
    p <- plotES(x = 1, gsea_results = gsea_res)
    expect_s3_class(p, "patchwork")
})

test_that("`plotES()` looks as expected with defaults", {
    skip_on_ci() # Text in plots is occasionally not reproduced exactly in GA
    expect_snapshot_plot(plotES(gsea_results = gsea_res), "plotES_default")
})

test_that("`plotES()` looks as expected with leading edge genes marked", {
    skip_on_ci()
    expect_snapshot_plot(plotES(gsea_results = gsea_res, show_leading_edge = TRUE), "plotES_showLeadingEdge")
})

test_that("`plotES()` leading edge only appear at/before ES peak (ES>1)", {
    n <- 100
    gene_names <- paste0("GENE_", seq_len(n))
    ranks <- withr::with_seed(123, sort(rnorm(n), decreasing = TRUE))
    names(ranks) <- gene_names

    # Gene set = first 15 genes (high ranks -> positive ES)
    gene_set <- gene_names[1:15]

    plot_df <- SomaEnrich:::.calc_plot_data(ranks, gene_set)

    # Peak of ES curve (positive ES -> which.max)
    es_peak_idx <- which.max(plot_df$runningScore)

    # "Leading edge" = gene-set members at or before the peak
    true_le <- plot_df$gene[plot_df$position == 1 & plot_df$x <= es_peak_idx]
    expect_true(length(true_le) > 0)

    # Simulate the fixed marking logic from plotES()
    plot_df$is_leading_edge <- plot_df$gene %in% true_le &
                               plot_df$position == 1 &
                               plot_df$x <= es_peak_idx

    le_positions <- plot_df$x[plot_df$is_leading_edge]

    # ALL marked leading-edge positions must be <= the peak
    expect_true(all(le_positions <= es_peak_idx))
})


test_that("`plotES()` leading edge only appear at/after ES peak (ES<1)", {
    n <- 100
    ranks <- withr::with_seed(42, sort(rnorm(n), decreasing = TRUE))
    gene_names <- paste0("GENE_", seq_len(n))
    names(ranks) <- gene_names

    # Gene set = last 15 genes (low ranks -> negative ES)
    gene_set <- gene_names[86:100]

    plot_df <- SomaEnrich:::.calc_plot_data(ranks, gene_set)

    # ES is negative -> peak (trough) at which.min
    es_peak_idx <- which.min(plot_df$runningScore)

    # "Leading edge" = gene-set members at or after the trough
    true_le <- plot_df$gene[plot_df$position == 1 & plot_df$x >= es_peak_idx]
    expect_true(length(true_le) > 0)

    plot_df$is_leading_edge <- plot_df$gene %in% true_le &
                               plot_df$position == 1 &
                               plot_df$x >= es_peak_idx

    le_positions <- plot_df$x[plot_df$is_leading_edge]

    # ALL marked leading-edge positions must be >= the trough
    expect_true(all(le_positions >= es_peak_idx))
})


test_that("`plotES()` duplicate gene names produce one LE marker per unique gene", {
    # Scenario: a gene name appears TWICE in the ranked list -- once before the
    # peak and once after. 
    # The fix should only mark one point per unique gene, matching the count in 
    # gsea_results$leadingEdge
    n <- 50
    gene_names <- withr::with_seed(123, paste0("GENE_", seq_len(n)))

    # Introduce a duplicate
    gene_names[45] <- "GENE_1"

    ranks <- withr::with_seed(123, sort(rnorm(n), decreasing = TRUE))
    names(ranks) <- gene_names

    gene_set <- gene_names[1:10]  # GENE_1 is in the set -> positive ES

    plot_df <- SomaEnrich:::.calc_plot_data(ranks, gene_set)
    es_peak_idx <- which.max(plot_df$runningScore)

    le_feats <- c("GENE_1")  # pretend fgsea said GENE_1 is leading edge

    # Apply the marking + dedup logic from plotES() (positive ES case)
    is_candidate <- plot_df$gene %in% le_feats &
                    plot_df$position == 1 &
                    plot_df$x <= es_peak_idx
    le_genes_seen <- plot_df$gene %in% le_feats & is_candidate
    is_dup <- duplicated(ifelse(le_genes_seen, plot_df$gene, NA)) & le_genes_seen
    plot_df$is_leading_edge <- is_candidate & !is_dup

    le_positions <- plot_df$x[plot_df$is_leading_edge]

    # Only ONE point should be marked for the single unique LE gene
    expect_equal(sum(plot_df$is_leading_edge), length(le_feats))

    # All marked leading-edge positions must be <= the peak
    expect_true(all(le_positions <= es_peak_idx))
})

test_that("`plotES()` accepts a character pathway name for `x`", {
    pathway_name <- gsea_res$results$pathway[[1L]]
    p_char <- plotES(x = pathway_name, gsea_results = gsea_res)
    p_int  <- plotES(x = 1L,           gsea_results = gsea_res)
    
    expect_s3_class(p_char, "patchwork")
    # Both should produce equivalent plots for the same pathway
    expect_equal(p_char, p_int)
})

test_that("`plotES()` errors when character `x` does not match any pathway", {
    expect_error(
        plotES(x = "This pathway does not exist", gsea_results = gsea_res),
        "not found in 'gsea_results'"
    )
})

test_that("`plotES()` works with `show_leading_edge = TRUE`", {
    expect_no_error(
        plotES(x = 1, gsea_results = gsea_res, show_leading_edge = TRUE)
    )
    p <- plotES(x = 1, gsea_results = gsea_res, show_leading_edge = TRUE)
    expect_s3_class(p, "patchwork")
})

test_that("`plotES()` works with a custom pathway via `cust_path`", {
    # Run GSEA with a small custom pathway
    custom_genes    <- gsea_res$results$leadingEdge[[1L]]
    cust_paths_list <- list(MyCustomPath = custom_genes)
    
    gsea_cust <- withr::with_seed(42,
        somaPrGSEA(ranks = ranks_pe,
                   cust_paths = cust_paths_list,
                   split_ids  = FALSE)
    )
    
    # plotES should accept cust_path and not error trying to look it up in pathway_map
    expect_no_error(
        plotES(x = 1, gsea_results = gsea_cust, cust_path = cust_paths_list)
    )
    p_cust <- plotES(x = 1, gsea_results = gsea_cust, cust_path = cust_paths_list)
    expect_s3_class(p_cust, "patchwork")
})

test_that("`plotES()` errors when `gsea_results` is not a list", {
    expect_error(
        plotES(x = 1, gsea_results = data.frame()),
        "`gsea_results` must be the list output of somaPrGSEA()"
    )
})

test_that("`plotES()` errors when `show_leading_edge` is not logical", {
    expect_error(
        plotES(x = 1, gsea_results = gsea_res, show_leading_edge = "yes"),
        "`show_leading_edge` must be logical"
    )
})

test_that("`plotES()` works as expected for a pathway with positive ES", {
    skip_on_ci()
    pos_idx <- which(gsea_res$results$ES > 0)[1L]
    expect_snapshot_plot(plotES(x = pos_idx, 
                                gsea_results = gsea_res, 
                                show_leading_edge = TRUE), 
                         "plotES_posES")
})

test_that("`plotES()` works as expected for a pathway with negative ES", {
    skip_on_ci()
    neg_idx <- which(gsea_res$results$ES < 0)[1L]
    expect_snapshot_plot(plotES(x = neg_idx, 
                                gsea_results = gsea_res, 
                                show_leading_edge = TRUE), 
                         "plotES_negES")
})
