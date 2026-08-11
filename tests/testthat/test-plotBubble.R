# Setup ----
deg <- head(t_tests$EntrezGeneSymbol, 50)
uni <- t_tests$EntrezGeneSymbol
res <- withr::with_seed(101, somaORA(features = deg, universe = uni))


# Testing ----
test_that("`plotBubble()` returns the expected plot with defaults", {
    skip_on_ci() # Text in plots is occasionally not reproduced exactly in GA
    expect_snapshot_plot(plotBubble(res), "plotBubble_default")
})

test_that("`plotBubble()` returns the expected ggplot object", {
    p <- plotBubble(res)

    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("`plotBubble()` errors when `ora_results` is not a data.frame", {
    expect_error(
        plotBubble(list(pathway = "test")),
        "`ora_results` must be a data.frame"
    )
    expect_error(
        plotBubble("not_a_df"),
        "`ora_results` must be a data.frame"
    )
})

test_that("`plotBubble()` respects `n_pathways` argument", {
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

test_that("`plotBubble(path_labels = 'id')` uses pathway IDs on y-axis", {
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
