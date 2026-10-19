#' Enrichment Bubble Plot
#'
#' Creates a bubble plot from ORA or GSEA results. The x-axis
#' metric depends on the class of the results: fold enrichment 
#' is used for [somaORA()] results, and NES for [somaPrGSEA()] 
#' results.
#'
#' @details
#' In ORA-derived bubble plots, bubble size corresponds
#' to the percentage of features in the feature set/pathway that are provided as
#' "interesting" features to [somaORA()].
#'
#' In GSEA-derived bubble plots, bubble size corresponds to the percentage of
#' features in the feature set/pathway that fall in the leading edge subset
#' returned by [somaPrGSEA()].
#'
#' In both types of plots, one bubble is plotted per pathway, and bubbles
#' are colored by adjusted p-value.
#'
#' @param x An object of class `somaORA` from [somaORA()], or of class
#'   `somaPrGSEA` from [somaPrGSEA()]. Unclassed ORA `data.frame`s and GSEA
#'   lists (containing `results` and `final_ranks` elements) are also accepted via
#'   the default method.
#' @param n_pathways The number of feature sets/pathways to be included in the
#'   plot. Default is 15.
#' @param font_size Numeric. Base font size for the plot. Default is 9.
#' @param bubble_color_low Character. Color to use for low (i.e. small)
#'   p-values. Default is `#FF9600` (orange).
#' @param bubble_color_high Character. Color to use for high (i.e. larger)
#'   p-values. Default is `#5164F8` (blue).
#' @param path_labels Character. String indicating the identifier to be used
#'   for each feature set. Options include `"name"` (the full feature set
#'   name) or `"id"` (the set identifier/accession number). Default is
#'   `"name"`. This sets the y-axis labels for both ORA and GSEA plots.
#' @param title Character or `NULL`. Optional plot title. Default is `NULL`.
#' @param grey_nonsig Numeric or logical. GSEA only. Bubbles with an adjusted
#'   p-value at or above this cutoff are filled grey. `TRUE` is equivalent to
#'   `0.05`, and `FALSE` disables greying. Default is `FALSE`. Numeric values
#'   (`0.1`, etc.) will set the adj. p-value cutoff to that value.
#' @param ... Additional arguments passed to methods.
#' @returns A ggplot object. Rows are the top `n_pathways` sets by adjusted
#'   p-value. For ORA, pathways are displayed sorted by adjusted p-value, the
#'   x-axis is fold enrichment, and bubble size is `% Overlap` (percentage of
#'   the pathway corresponding to interesting features provided to
#'   [somaORA()]). For GSEA, pathways are sorted by adjusted p-value, the
#'   x-axis is NES, and bubble size is leading-edge percentage.
#' @author Amanda Hiser
#' @examples
#' # --- ORA bubble plot ---
#' deg <- head(t_tests$EntrezGeneSymbol, 50)
#' bg <- SomaDataIO::getAnalyteInfo(example_data_11k)$EntrezGeneSymbol
#' ora_res <- somaORA(features = deg,
#'                    universe = bg)
#' plotBubble(ora_res)
#'
#' # Display fewer feature sets
#' plotBubble(ora_res, n_pathways = 10)
#'
#' # Use GO or MSigDb accession identifiers, instead of full names
#' plotBubble(ora_res, path_labels = "id")
#'
#' # Add a title
#' plotBubble(ora_res, title = "Overrepresentation Analysis (ORA) Results")
#'
#' # Other plot customizations can be performed using ggplot2
#' plotBubble(ora_res, title = "ORA Results") +
#'    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
#'                   panel.border = ggplot2::element_rect(colour = "black",
#'                                               fill = NA, linewidth = 1))
#'
#' \dontrun{
#' # --- GSEA bubble plot ---
#' ranks <- prepareRanks(stats = t_tests$log2_fc,
#'                       features = t_tests$EntrezGeneSymbol,
#'                       verbose = FALSE)
#' gsea_res <- somaPrGSEA(ranks = ranks, verbose = FALSE)
#' plotBubble(gsea_res, title = "GSEA Results")
#'
#' # Grey out bubbles with adjusted p-value >= 0.1
#' plotBubble(gsea_res, grey_nonsig = 0.1)
#' }
#' @importFrom stringr str_to_sentence
#' @importFrom ggplot2 ggplot aes geom_point scale_size scale_fill_gradient
#' @importFrom ggplot2 labs theme element_text element_blank element_line
#' @importFrom ggplot2 guides guide_legend guide_colorbar scale_y_discrete
#' @importFrom scales label_wrap
#' @rdname plotBubble
#' @export
plotBubble <- function(x, ...) {
  UseMethod("plotBubble")
}


#' @rdname plotBubble
#' @export
plotBubble.somaORA <- function(x,
                               n_pathways = 15,
                               font_size = 9,
                               bubble_color_low = "#FF9600",
                               bubble_color_high = "#5164F8",
                               path_labels = c("name", "id"),
                               title = NULL,
                               ...) {

  path_labels <- match.arg(path_labels)

  # Calculate percentage of pathway covered by DE features
  plot_df <- x
  plot_df$pctDE <- plot_df$overlap / plot_df$size * 100

  # Sort by adj. p, pctDE, and pathway name (keeps order stable in case of ties)
  plot_df <- plot_df[order(plot_df$padj, -plot_df$pctDE, plot_df$pathway), ]

  # Restrict to number of results specified by the user
  n_pathways <- min(n_pathways, nrow(plot_df))
  plot_df <- plot_df[seq_len(n_pathways), ]

  # Create new ordered factor to use for plotting
  if ( path_labels == "name" ) {
    plot_df$category <- plot_df$pathway
  } else {
    plot_df$category <- plot_df$pathway_id
  }

  plot_df$category <- factor(plot_df$category,
                             levels = plot_df$category[order(-plot_df$padj,
                                                             plot_df$pctDE,
                                                             -seq_len(nrow(plot_df)))])

  yax_lab <- unique(group_code_lookup[group_code_lookup$group_code %in% plot_df$resource_code, "label"])

  p <- ggplot(plot_df, aes(x = foldEnrichment,
                           y = category,
                           size  = pctDE,
                           fill  = padj)) +
    geom_point(alpha = 0.75, shape = 21, color = "black") +
    scale_size(range = c(2, 12)) +
    scale_fill_gradient(limits = c(0, max(plot_df$padj)),
                        low = bubble_color_low,
                        high = bubble_color_high) +
    guides(fill = guide_colorbar(order = 1, reverse = TRUE),
           color = "none",
           size = guide_legend(order = 2,
                               override.aes = list(fill = "grey10"))) +
    labs(x     = "\nFold Enrichment",
         y     = paste(yax_lab, "\n"),
         size  = "% Overlap",
         fill  = "Adjusted\np-value",
         title = title) +
    theme(legend.position = "right",
          axis.text.y = element_text(size = font_size, angle = 0, hjust = 1),
          axis.title  = element_text(size = font_size + 2),
          axis.ticks  = element_blank(),
          panel.background  = element_blank(),
          plot.background   = element_blank(),
          panel.border      = element_blank(),
          panel.grid.major  = element_line(colour = "grey92"),
          panel.grid.minor  = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_y_discrete(labels = function(x) {
      scales::label_wrap(30)(
        if ( path_labels == "name" ) stringr::str_to_sentence(x) else x
      )
    })

  return(p)
}


#' @rdname plotBubble
#' @export
plotBubble.somaPrGSEA <- function(x,
                                  n_pathways = 15,
                                  font_size = 9,
                                  bubble_color_low = "#FF9600",
                                  bubble_color_high = "#5164F8",
                                  path_labels = c("name", "id"),
                                  title = NULL,
                                  grey_nonsig = FALSE,
                                  ...) {

  path_labels <- match.arg(path_labels)

  if ( isTRUE(grey_nonsig) ) {
    grey_nonsig <- 0.05
  } else if ( isFALSE(grey_nonsig) ) {
    grey_nonsig <- Inf
  } else if ( !is.numeric(grey_nonsig) || length(grey_nonsig) != 1L ) {
    stop("`grey_nonsig` must be TRUE, FALSE, or a single number.", call. = FALSE)
  }

  plot_df <- x$results
  plot_df$pct_le <- lengths(plot_df$leadingEdge) / plot_df$final_set_size * 100

  # Select top pathways by adj. p, magnitude of enrichment, and set size
  plot_df <- plot_df[order(plot_df$padj, -abs(plot_df$NES),
                           plot_df$final_set_size), ]

  n_pathways <- min(n_pathways, nrow(plot_df))
  plot_df <- plot_df[seq_len(n_pathways), ]

  if ( path_labels == "name" ) {
    plot_df$category <- plot_df$pathway
  } else {
    plot_df$category <- plot_df$pathway_id
  }

  # Most significant pathways at top. First factor level is drawn at the bottom.
  plot_df$category <- factor(plot_df$category,
                             levels = plot_df$category[order(-plot_df$padj,
                                                             plot_df$pct_le,
                                                             -seq_len(nrow(plot_df)))])

  plot_df$fill_padj <- replace(plot_df$padj, plot_df$padj >= grey_nonsig, NA)

  # Labels created using internal group_code_lookup obj (in R/sysdata.rda)
  yax_lab <- unique(
    group_code_lookup[group_code_lookup$group_code %in% plot_df$resource_code, "label"]
  )

  p <- ggplot(plot_df, aes(x = NES,
                           y = category,
                           size = pct_le,
                           fill = fill_padj)) +
    geom_point(alpha = 0.75, shape = 21, color = "black") +
    scale_size(range = c(2, 12)) +
    scale_fill_gradient(limits = c(0, NA),
                        low = bubble_color_low,
                        high = bubble_color_high,
                        na.value = "#ADADAD") +
    guides(fill = guide_colorbar(order = 1, reverse = TRUE),
           color = "none",
           size = guide_legend(order = 2,
                               override.aes = list(fill = "grey10"))) +
    labs(x     = "\nNormalized Enrichment Score",
         y     = paste(yax_lab, "\n"),
         size  = "% LE Genes in Pathway",
         fill  = "Adjusted\np-value",
         title = title) +
    theme(legend.position = "right",
          axis.text.y = element_text(size = font_size, angle = 0, hjust = 1),
          axis.title  = element_text(size = font_size + 2),
          axis.ticks  = element_blank(),
          panel.background = element_blank(),
          plot.background  = element_blank(),
          panel.border     = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_y_discrete(labels = function(x) {
      scales::label_wrap(30)(
        if ( path_labels == "name" ) stringr::str_to_sentence(x) else x
      )
    })

  return(p)
}


#' @rdname plotBubble
#' @export
plotBubble.default <- function(x, ...) {

  if ( is.data.frame(x) &&
       all(c("foldEnrichment", "overlap", "size", "padj", "pathway") %in% names(x)) ) {
    class(x) <- c("somaORA", setdiff(class(x), "somaORA"))
    return(plotBubble(x, ...))
  }

  if ( is.list(x) &&
       all(c("results", "final_ranks") %in% names(x)) &&
       is.data.frame(x$results) &&
       all(c("NES", "padj", "leadingEdge", "final_set_size") %in% names(x$results)) ) {
    class(x) <- c("somaPrGSEA", setdiff(class(x), "somaPrGSEA"))
    return(plotBubble(x, ...))
  }

  stop(
    "`x` must be the output of somaORA() or somaPrGSEA().",
    call. = FALSE
  )
}
