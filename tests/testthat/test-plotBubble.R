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
test_that("`plotBubble()` returns the expected plot with defaults", {
    skip_on_ci() # Text in plots is occasionally not reproduced exactly in GA
    expect_snapshot_plot(plotBubble(res), "plotBubble_default")
})

test_that("`plotBubble()` returns the expected ggplot object for ORA", {
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

    # Extracting actual number of rows plotted
    n5  <- nrow(ggplot2::ggplot_build(p5)$data[[1]])
    n10 <- nrow(ggplot2::ggplot_build(p10)$data[[1]])

    expect_equal(n5, 5L)
    expect_equal(n10, 10L)
})

test_that("`plotBubble()` does not error when `n_pathways` exceeds nrow(results)", {
    # Previously would produce NA rows in plot_df -> ggplot error
    n_large <- nrow(res) + 100L
    expect_no_error(plotBubble(res, n_pathways = n_large))

    # Should plot all available rows, not NA-padded ones
    p <- plotBubble(res, n_pathways = n_large)
    n_points <- nrow(ggplot2::ggplot_build(p)$data[[1]])
    expect_equal(n_points, nrow(res))
})

test_that("`plotBubble(path_labels = 'id')` uses pathway IDs on y-axis for ORA", {
    p_name <- plotBubble(res, path_labels = "name")
    p_id   <- plotBubble(res, path_labels = "id")

    # These should differ from each other
    expect_false(identical(p_name, p_id))

    # Build the id plot and confirm y-axis labels look like IDs, not names
    built_id <- ggplot2::ggplot_build(p_id)
    y_labels  <- built_id$layout$panel_params[[1]]$y$get_labels()
    y_labels  <- y_labels[!is.na(y_labels)]
    expect_true(any(grepl("^GO:", y_labels)))
})

test_that("`plotBubble()` `path_labels` defaults to 'name' and errors on bad value", {
    expect_no_error(plotBubble(res))  # Default "name"
    expect_error(plotBubble(res, path_labels = "invalid"),
                 "'arg' should be one of")
})

test_that("`plotBubble()` custom bubble colors are accepted without error", {
    expect_no_error(
        plotBubble(res, bubble_color_low = "red", bubble_color_high = "green")
    )
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

test_that("`plotBubble()` returns a ggplot object for GSEA", {
    p <- plotBubble(gsea_res)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("`plotBubble()` GSEA axis and legend labels are as expected", {
    labs <- plotBubble(gsea_res)$labels

    expect_equal(labs$x, "\nNES")
    expect_equal(labs$y, "GO Biological Process \n")
    expect_equal(labs$size, "Leading\nedge %")
    expect_equal(labs$fill, "Adjusted\np-value")
})

test_that("`plotBubble()` respects `n_pathways` argument for GSEA", {
    p5  <- plotBubble(gsea_res, n_pathways = 5)
    p10 <- plotBubble(gsea_res, n_pathways = 10)

    n5  <- nrow(ggplot2::ggplot_build(p5)$data[[1]])
    n10 <- nrow(ggplot2::ggplot_build(p10)$data[[1]])

    expect_equal(n5, 5L)
    expect_equal(n10, 10L)
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
    expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 5L)
})

test_that("`plotBubble()` calculates leading edge % from `leadingEdge` and `final_set_size`", {
    p <- plotBubble(gsea_res)
    expected <- lengths(p$data$leadingEdge) / p$data$final_set_size * 100

    expect_equal(p$data$pct_le, expected)
    expect_true(all(p$data$pct_le > 0 & p$data$pct_le <= 100))
})

test_that("`plotBubble()` facets GSEA pathways by NES direction", {
    p <- plotBubble(gsea_res)
    built <- ggplot2::ggplot_build(p)
    layout <- built$layout$layout

    expect_s3_class(p$facet, "FacetGrid")
    expect_setequal(as.character(layout$direction), c("Negative NES", "Positive NES"))

    pts <- built$data[[1]]
    panel_dir <- as.character(layout$direction[match(pts$PANEL, layout$PANEL)])
    expect_true(all(pts$x[panel_dir == "Positive NES"] > 0))
    expect_true(all(pts$x[panel_dir == "Negative NES"] < 0))
})

test_that("`plotBubble()` does not facet GSEA plots with a single NES direction", {
    for ( p in list(plotBubble(gsea_pos), plotBubble(gsea_neg)) ) {
        expect_s3_class(p$facet, "FacetNull")
        expect_s3_class(p$theme$panel.border, "element_blank")
    }
})

test_that("`plotBubble()` does not facet GSEA plots when `n_pathways` yields one NES direction", {
    p <- plotBubble(gsea_res, n_pathways = 1)
    expect_s3_class(p$facet, "FacetNull")
})

test_that("`plotBubble()` sorts GSEA pathways by adj. p within each NES direction", {
    p <- plotBubble(gsea_res)
    df <- p$data
    # Top of plot is the last factor level
    df <- df[order(df$category, decreasing = TRUE), ]

    expect_true(all(diff(as.integer(df$direction)) <= 0))
    for ( dir in levels(df$direction) ) {
        expect_false(is.unsorted(df$padj[df$direction == dir]))
    }
})

test_that("`plotBubble(path_labels = 'id')` uses pathway IDs on y-axis for GSEA", {
    p_name <- plotBubble(gsea_res, n_pathways = 5, path_labels = "name")
    p_id   <- plotBubble(gsea_res, n_pathways = 5, path_labels = "id")

    expect_false(identical(p_name, p_id))

    built_id <- ggplot2::ggplot_build(p_id)
    y_labels <- built_id$layout$panel_params[[1]]$y$get_labels()
    y_labels <- y_labels[!is.na(y_labels)]
    expect_true(any(grepl("^(GO:|M[0-9])", y_labels)))
})

test_that("`plotBubble()` `path_labels` errors on bad value for GSEA", {
    expect_error(plotBubble(gsea_res, path_labels = "invalid"),
                 "'arg' should be one of")
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
