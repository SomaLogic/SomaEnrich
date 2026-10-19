# Setup ----
deg <- head(t_tests$EntrezGeneSymbol, 50)
uni <- t_tests$EntrezGeneSymbol
res <- withr::with_seed(101, somaORA(features = deg, universe = uni))

ranks <- t_tests$t_stat
names(ranks) <- t_tests$EntrezGeneSymbol
ranks <- sort(ranks, decreasing = TRUE)
gsea_res <- withr::with_seed(101, suppressWarnings(somaPrGSEA(ranks = ranks)))

gsea_pos <- gsea_res
gsea_pos$results <- gsea_pos$results[gsea_pos$results$NES > 0, ]

gsea_neg <- gsea_res
gsea_neg$results <- gsea_neg$results[gsea_neg$results$NES < 0, ]


# Testing ----
# ORA method ----
test_that("`plotBubble()` returns the expected ORA plot with defaults", {
    skip_on_ci() # Text in plots is occasionally not reproduced exactly in GA
    expect_snapshot_plot(plotBubble(res), "plotBubble_default")
})

test_that("`plotBubble(path_labels = 'id')` returns the expected ORA plot", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(res, n_pathways = 10, path_labels = "id"),
                         "plotBubble_ora_ids")
})

test_that("`plotBubble()` returns the expected ORA plot with custom styling", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(res,
                                    n_pathways = 10,
                                    font_size = 12,
                                    bubble_color_low = "red",
                                    bubble_color_high = "grey50"),
                         "plotBubble_ora_custom")
})

test_that("`plotBubble()` returns the expected ORA plot with a title", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(res, n_pathways = 10, title = "ORA Bubble Plot"),
                         "plotBubble_ora_title")
})

test_that("`plotBubble()` returns a ggplot object for ORA", {
    p <- plotBubble(res)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("`plotBubble()` errors when `x` is not ORA or GSEA output", {
    expect_error(
        plotBubble(list(pathway = "test")),
        "`x` must be the output of somaORA() or somaPrGSEA()",
        fixed = TRUE
    )
    expect_error(
        plotBubble("not_a_df"),
        "`x` must be the output of somaORA() or somaPrGSEA()",
        fixed = TRUE
    )
})

test_that("`plotBubble()` respects `n_pathways` argument for ORA", {
    p5  <- plotBubble(res, n_pathways = 5)
    p10 <- plotBubble(res, n_pathways = 10)

    expect_equal(nrow(p5$data), 5L)
    expect_equal(nrow(p10$data), 10L)
})

test_that("`plotBubble()` does not error when `n_pathways` exceeds nrow(results)", {
    # Previously would produce NA rows in plot_df -> ggplot error
    n_large <- nrow(res) + 100L
    expect_no_error(plotBubble(res, n_pathways = n_large))

    p <- plotBubble(res, n_pathways = n_large)
    expect_equal(nrow(p$data), nrow(res))
})

test_that("`plotBubble()` `path_labels` defaults to 'name' and errors on bad value", {
    expect_no_error(plotBubble(res))  # Default "name"
    expect_error(plotBubble(res, path_labels = "invalid"),
                 "'arg' should be one of")
})

# GSEA method ----
test_that("`plotBubble()` returns the expected GSEA plot with defaults", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_res), "plotBubble_gsea_default")
})

test_that("`plotBubble()` returns the expected GSEA plot with only positive NES", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_pos), "plotBubble_gsea_posNES")
})

test_that("`plotBubble()` returns the expected GSEA plot with only negative NES", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_neg), "plotBubble_gsea_negNES")
})

test_that("`plotBubble(path_labels = 'id')` returns the expected GSEA plot", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_res, n_pathways = 10, path_labels = "id"),
                         "plotBubble_gsea_ids")
})

test_that("`plotBubble()` returns the expected GSEA plot with custom styling", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_res,
                                    n_pathways = 15,
                                    font_size = 12,
                                    bubble_color_low = "red",
                                    bubble_color_high = "grey50"),
                         "plotBubble_gsea_custom")
})

test_that("`plotBubble()` returns the expected GSEA plot with `grey_nonsig`", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_res, n_pathways = 15, grey_nonsig = 0.05),
                         "plotBubble_gsea_grey_nonsig")
})

test_that("`plotBubble()` returns the expected GSEA plot with a title", {
    skip_on_ci()
    expect_snapshot_plot(plotBubble(gsea_res, n_pathways = 10, title = "GSEA Bubble Plot"),
                         "plotBubble_gsea_title")
})

test_that("`plotBubble()` returns a ggplot object for GSEA", {
    p <- plotBubble(gsea_res)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("`plotBubble()` respects `n_pathways` argument for GSEA", {
    p5  <- plotBubble(gsea_res, n_pathways = 5)
    p10 <- plotBubble(gsea_res, n_pathways = 10)

    expect_equal(nrow(p5$data), 5L)
    expect_equal(nrow(p10$data), 10L)
})

test_that("`plotBubble()` selects the top GSEA pathways by adjusted p-value", {
    p <- plotBubble(gsea_res, n_pathways = 10)
    top_padj <- sort(gsea_res$results$padj)[seq_len(10)]

    expect_equal(sort(p$data$padj), top_padj)
})

test_that("`plotBubble()` does not error when `n_pathways` exceeds nrow(results) for GSEA", {
    small <- gsea_res
    small$results <- small$results[order(small$results$padj), ][seq_len(5), ]

    expect_no_error(plotBubble(small, n_pathways = 100))
    p <- plotBubble(small, n_pathways = 100)
    expect_equal(nrow(p$data), 5L)
})

test_that("`plotBubble()` calculates leading edge % from `leadingEdge` and `final_set_size`", {
    p <- plotBubble(gsea_res)
    expected <- lengths(p$data$leadingEdge) / p$data$final_set_size * 100

    expect_equal(p$data$pct_le, expected)
    expect_true(all(p$data$pct_le > 0 & p$data$pct_le <= 100))
})

test_that("`plotBubble()` sorts GSEA pathways by adjusted p-value", {
    p <- plotBubble(gsea_res, n_pathways = 10)
    # Top of plot is the last factor level (most significant / lowest padj)
    df <- p$data[order(p$data$category, decreasing = TRUE), ]
    expect_false(is.unsorted(df$padj))
})

test_that("`plotBubble()` `path_labels` errors on bad value for GSEA", {
    expect_error(plotBubble(gsea_res, path_labels = "invalid"),
                 "'arg' should be one of")
})

test_that("`plotBubble()` `grey_nonsig` treats TRUE as 0.05 and FALSE as no cutoff", {
    expect_equal(plotBubble(gsea_res, grey_nonsig = TRUE)$data,
                 plotBubble(gsea_res, grey_nonsig = 0.05)$data)
    p_off <- plotBubble(gsea_res, grey_nonsig = FALSE)
    expect_equal(p_off$data$fill_padj, p_off$data$padj)
})

test_that("`plotBubble()` marks GSEA bubbles at or above `grey_nonsig` as NA fill", {
    p <- plotBubble(gsea_res, grey_nonsig = 0.001)
    nonsig <- p$data$padj >= 0.001

    expect_true(any(nonsig) && !all(nonsig))
    expect_true(all(is.na(p$data$fill_padj[nonsig])))
    expect_equal(p$data$fill_padj[!nonsig], p$data$padj[!nonsig])
})

test_that("`plotBubble()` errors on invalid `grey_nonsig` values", {
    expect_error(plotBubble(gsea_res, grey_nonsig = "yes"),
                 regexp = "`grey_nonsig` must be TRUE, FALSE, or a single number.",
                 fixed = TRUE)
})

# Default method ----
test_that("`plotBubble.default()` accepts an unclassed ORA data.frame", {
    bare <- res
    class(bare) <- "data.frame"
    p <- plotBubble(bare)

    expect_s3_class(p, "ggplot")
    expect_equal(p$data, plotBubble(res)$data)
})

test_that("`plotBubble.default()` accepts an unclassed GSEA list", {
    bare <- unclass(gsea_res)
    p <- plotBubble(bare)

    expect_s3_class(p, "ggplot")
    expect_equal(p$data, plotBubble(gsea_res)$data)
})

test_that("`plotBubble.default()` passes arguments through to the dispatched method", {
    bare <- unclass(gsea_res)
    p <- plotBubble(bare, n_pathways = 5, path_labels = "id")

    expect_equal(p$data, plotBubble(gsea_res, n_pathways = 5, path_labels = "id")$data)
})

test_that("`plotBubble.default()` errors when GSEA-like input is missing required columns", {
    bad_gsea <- unclass(gsea_res)
    bad_gsea$results$leadingEdge <- NULL
    expect_error(plotBubble(bad_gsea),
                 "`x` must be the output of somaORA() or somaPrGSEA().",
                 fixed = TRUE)

    bad_ora <- res
    class(bad_ora) <- "data.frame"
    bad_ora$foldEnrichment <- NULL
    expect_error(plotBubble(bad_ora),
                 "`x` must be the output of somaORA() or somaPrGSEA().",
                 fixed = TRUE)
})
