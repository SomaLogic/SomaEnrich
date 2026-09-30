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
#' returned by [somaPrGSEA()]. GSEA plots are split into side-by-side panels
#' for negative and positive NES, each with its own x-axis range. If all
#' plotted pathways share the same NES direction, a single unfaceted panel is
#' shown instead.
#'
#' In both types of plots, one bubble is plotted per pathway, and bubbles
#' are colored by adjusted p-value.
#'
#' @param x An object of class `somaORA` from [somaORA()], or of class
#'   `somaPrGSEA` from [somaPrGSEA()]. Unclassed ORA `data.frame`s and GSEA
#'   lists (containing `results` and `final_ranks` elements) are also accepted via
#'   the default method.
#' @param n_pathways The number of feature sets/pathways to be included in the
#'   plot. Default is 25.
#' @param font_size Numeric. Base font size for the plot. Default is 9.
#' @param bubble_color_low Character. Color to use for low (i.e. small)
#'   p-values. Default is `"orange"`.
#' @param bubble_color_high Character. Color to use for high (i.e. larger)
#'   p-values. Default is `"blue"`.
#' @param path_labels Character. String indicating the identifier to be used
#'   for each feature set. Options include `"name"` (the full feature set
#'   name) or `"id"` (the set identifier/accession number). Default is
#'   `"name"`. This sets the y-axis labels for both ORA and GSEA plots.
#' @param ... Additional arguments passed to methods.
#' @returns A ggplot object. Rows are the top `n_pathways` sets by adjusted
#'   p-value. For ORA, pathways are displayed sorted by adjusted p-value, the
#'   x-axis is fold enrichment, and bubble size is `% Overlap` (percentage of
#'   the pathway corresponding to interesting features provided to
#'   [somaORA()]). For GSEA, pathways with positive NES are displayed above
#'   those with negative NES, each group sorted by adjusted p-value; the
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
#' # Other plot customizations can be performed using ggplot2
#' plotBubble(ora_res) +
#'    ggplot2::ggtitle("Overrepresentation Analysis (ORA) Results") +
#'    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
#'                   panel.border = ggplot2::element_rect(colour = "black", 
#'                                               fill = NA, linewidth = 1))
#'
#' \dontrun{
#' # --- GSEA bubble plot ---
#' ranks <- prepareRanks(stats = t_tests$t_stat,
#'                       features = t_tests$EntrezGeneSymbol)
#' gsea_res <- somaPrGSEA(ranks = ranks)
#' plotBubble(gsea_res)
#' # Only one panel will be generated if all NES values are positive or negative
#' plotBubble(gsea_res, n_pathways = 10)
#' }
#' @importFrom stringr str_to_sentence
#' @importFrom ggplot2 ggplot aes geom_point scale_size scale_fill_gradient
#' @importFrom ggplot2 labs theme element_text element_blank element_line
#' @importFrom ggplot2 guides guide_legend guide_colorbar scale_y_discrete
#' @importFrom ggplot2 facet_grid vars scale_x_continuous expansion element_rect
#' @importFrom scales label_wrap breaks_extended
#' @importFrom rlang .data
#' @rdname plotBubble
#' @export
plotBubble <- function(x, ...) {
  UseMethod("plotBubble")
}


#' @rdname plotBubble
#' @export
plotBubble.somaORA <- function(x,
                               n_pathways = 25,
                               font_size = 9,
                               bubble_color_low = "orange",
                               bubble_color_high = "blue",
                               path_labels = c("name", "id"),
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
    plot_df$category <- stringr::str_to_sentence(plot_df$pathway)
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
         fill = "Adjusted\np-value") +
    theme(legend.position = "right",
          axis.text.y = element_text(size = font_size, angle = 0, hjust = 1),
          axis.title  = element_text(size = font_size + 2),
          axis.ticks  = element_blank(),
          panel.background  = element_blank(),
          plot.background   = element_blank(),
          panel.border      = element_blank(),
          panel.grid.major  = element_line(colour = "grey92"),
          panel.grid.minor  = element_blank()) +
    scale_y_discrete(labels = scales::label_wrap(30))

  return(p)
}


#' @rdname plotBubble
#' @export
plotBubble.somaPrGSEA <- function(x,
                                  n_pathways = 25,
                                  font_size = 9,
                                  bubble_color_low = "orange",
                                  bubble_color_high = "blue",
                                  path_labels = c("name", "id"),
                                  ...) {

  path_labels <- match.arg(path_labels)

  plot_df <- x$results
  plot_df$pct_le <- lengths(plot_df$leadingEdge) / plot_df$final_set_size * 100

  # Select top pathways by adj. p, then |NES| (stable ties by pathway name)
  plot_df <- plot_df[order(plot_df$padj, -abs(plot_df$NES), plot_df$pathway), ]

  n_pathways <- min(n_pathways, nrow(plot_df))
  plot_df <- plot_df[seq_len(n_pathways), ]

  if ( path_labels == "name" ) {
    plot_df$category <- stringr::str_to_sentence(plot_df$pathway)
  } else {
    plot_df$category <- plot_df$pathway_id
  }

  plot_df$direction <- factor(ifelse(plot_df$NES > 0, "Positive NES", "Negative NES"),
                              levels = c("Negative NES", "Positive NES"))

  # Positive NES pathways above negative, each block sorted by adj. p (most
  # significant at top). First factor level is drawn at the bottom.
  plot_df$category <- factor(plot_df$category,
                             levels = plot_df$category[order(plot_df$direction,
                                                             -plot_df$padj,
                                                             abs(plot_df$NES),
                                                             -seq_len(nrow(plot_df)))])

  yax_lab <- unique(group_code_lookup[group_code_lookup$group_code %in% plot_df$resource_code, "label"])

  p <- ggplot(plot_df, aes(x = NES,
                           y = category,
                           size = pct_le,
                           fill = padj)) +
    geom_point(alpha = 0.75, shape = 21, color = "black") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4),
                       expand = expansion(mult = 0.2)) +
    scale_size(range = c(2, 12)) +
    scale_fill_gradient(limits = c(0, max(plot_df$padj)),
                        low = bubble_color_low,
                        high = bubble_color_high) +
    guides(fill = guide_colorbar(order = 1, reverse = TRUE),
           color = "none",
           size = guide_legend(order = 2,
                               override.aes = list(fill = "grey10"))) +
    labs(x = "\nNES",
         y = paste(yax_lab, "\n"),
         size = "Leading\nedge %",
         fill = "Adjusted\np-value") +
    theme(legend.position = "right",
          axis.text.y = element_text(size = font_size, angle = 0, hjust = 1),
          axis.title = element_text(size = font_size + 2),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.background  = element_blank(),
          panel.border     = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_blank()) +
    scale_y_discrete(labels = scales::label_wrap(30))

  if ( nlevels(droplevels(plot_df$direction)) > 1L ) {
    p <- p +
      facet_grid(cols = vars(direction), scales = "free_x") +
      theme(panel.border     = element_rect(colour = "grey92", fill = NA),
            strip.background = element_rect(fill = "grey92", colour = NA),
            strip.text       = element_text(size = font_size + 1))
  }

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
