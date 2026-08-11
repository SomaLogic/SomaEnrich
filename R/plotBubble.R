#' ORA Bubble Plot
#' 
#' Creates a bubble plot to display results of ORA. Bubble size corresponds
#' to the percentage of features in the feature set/pathway that are provided as 
#' "interesting" features to [somaORA()]. Bubbles are colored by adjusted 
#' p-value. One bubble is plotted per pathway.
#' 
#' @inheritParams plotES
#' @param ora_results A `data.frame` of ORA results from [somaORA()].
#' @param n_pathways The number of feature sets/pathways to be included in the 
#'   plot. Default is 25.
#' @param bubble_color_low Character. Color to use for low (i.e. small) 
#'   p-values. Default is "orange".
#' @param bubble_color_high Character. Color to use for high (i.e. larger) 
#'   p-values. Default is "blue".
#' @param path_labels Character. String indicating the identifier to be used 
#'   for each feature set (on the y axis). Options include "name" 
#'   (the full feature set name) or "id" (the set identifier/accession number). 
#'   Default is "name".
#' @returns A ggplot object of the ORA bubble plot. Rows are sorted by
#'   the adjusted p-value and "% Overlap", i.e. percentage of the pathway 
#'   corresponding to "interesting" features (those provided to the `features` 
#'   argument of [somaORA()]).
#' @author Amanda Hiser
#' @examples
#' deg <- head(t_tests$EntrezGeneSymbol, 50)
#' bg <- SomaDataIO::getAnalyteInfo(example_data_11k)$EntrezGeneSymbol
#' res <- somaORA(features = deg, 
#'                universe = bg)
#' plotBubble(res)
#' 
#' # Display fewer feature sets
#' plotBubble(res, n_pathways = 10)
#' 
#' # Use GO or MSigDb accession identifiers, instead of full names
#' plotBubble(res, path_labels = "id")
#'    
#' # Other plot customizations can be performed using ggplot2
#' plotBubble(res) +
#'    ggplot2::ggtitle("Overrepresentation Analysis (ORA) Results") +
#'    ggplot2::theme(legend.position = "none")
#' @importFrom stringr str_to_sentence
#' @importFrom ggplot2 ggplot aes geom_point scale_size scale_fill_gradient
#' @importFrom ggplot2 labs theme element_text element_blank
#' @importFrom ggplot2 guides guide_legend guide_colorbar scale_y_discrete
#' @importFrom scales label_wrap
#' @importFrom rlang .data
#' @export
plotBubble <- function(ora_results, 
                       n_pathways = 25,
                       font_size = 9,
                       bubble_color_low = "orange",
                       bubble_color_high = "blue",
                       path_labels = c("name", "id")) {
    
    path_labels <- match.arg(path_labels)
    
    stopifnot("`ora_results` must be a data.frame" = is.data.frame(ora_results))

    # Calculate percentage of pathway covered by DE features
    plot_df <- ora_results
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
    
    ggplot(plot_df, aes(x = foldEnrichment, 
                        y = category, 
                        size  = pctDE,
                        fill  = padj)) +
      geom_point(alpha = 0.75, shape = 21, color = "black") +
      scale_size(range = c(2, 12)) +
      scale_fill_gradient(limits = c(0, max(plot_df$padj)),
                                   low = bubble_color_low,
                                   high = bubble_color_high) +
      guides(fill = guide_colorbar(order = 1, reverse = TRUE), # Continuous gradient bar for fill
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
}
