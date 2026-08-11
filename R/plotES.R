#' GSEA Enrichment Plot with Leading Edge Annotation
#'
#' Creates a GSEA enrichment plot, with the option to annotate leading edge 
#' features on the running enrichment score curve.
#'
#' @param x Character or integer. The pathway name _or_ row index of 
#'   `gsea_results` to use for creating the plot and extracting leading edge 
#'   features. Default is 1 (first row). If a character, must match a pathway
#'   name in the "pathway_name" column of `gsea_results`.
#' @param gsea_results A `data.frame` of GSEA results from [somaPrGSEA()]. 
#'   The `leadingEdge` column will be used to identify and annotate leading 
#'   edge features on the plot. The `resource_id` column will be used to
#'   retrieve pathway members. If a custom pathway was used for GSEA, the
#'   pathway must be provided to `cust_path`.
#' @param cust_path Optional. A vector of gene identifiers (gene symbols 
#'    or Entrez gene IDs) or SomaScan `AptNames` representing custom gene sets 
#'    or functional groups of interest. Should be provided if the `cust_paths`
#'    argument was used to generate results via [somaPrGSEA()].
#' @param show_leading_edge Logical. Should leading edge features be highlighted
#'   on the plot? If TRUE, leading edge features are extracted from 
#'   `gsea_results` and annotated with colored dots. Default is TRUE.
#' @param show_path_stats Logical. Should GSEA result statistics (ES, NES, 
#'   and path size) be annotated on the plot? Default is TRUE.
#' @param stats_color Character. Color of the stats annotation from
#'   `show_path_stats`.
#' @param enrichment_score_color Character. Color of the running enrichment 
#'   score line. Default is "green".
#' @param font_size Numeric. Base font size for the plot. Default is 12.
#' @param font_size Numeric. Base font size for the plot. Default is 12.
#' @param leading_edge_color Character. Color of the leading edge gene dots.
#'   Default is "red".
#' @param leading_edge_size Numeric. Size of the leading edge gene dots.
#'   Default is 2.
#' @param leading_edge_alpha Numeric. Transparency of the leading edge dots. 
#'   Default is 0.8.
#' @returns A ggplot object (combined via [patchwork]) displaying the 
#'   GSEA enrichment plot with optional leading edge gene annotation.
#' @details
#'   This function creates a plot with three panels:
#'   \enumerate{
#'     \item \strong{Enrichment Score Panel}: Shows the running enrichment score
#'       with leading edge genes optionally highlighted as colored dots.
#'     \item \strong{Gene Hits Panel}: Shows vertical lines indicating where
#'       genes in the gene set appear in the ranked list.
#'     \item \strong{Ranked List Metric Panel}: Shows the distribution of the
#'       ranking metric across all genes.
#'  }
#'   The leading edge genes are the core subset that contribute most to the
#'   enrichment signal, appearing at or before the point where the running
#'   enrichment score reaches its maximum deviation from zero. For a positive 
#'   ES, the leading edge subset is the set of members that appear in the 
#'   ranked list _prior_ to the peak score. For a negative ES, the leading edge 
#'   subset will appear subsequent to the peak score.
#' @references
#'   Subramanian, Tamayo, et al. Gene set enrichment analysis: A knowledge-based
#'   approach for interpreting genome-wide expression profiles. 
#'   Proc Natl Acad Sci USA. 102(43):15545-50 (2005). 
#'   https://doi.org/10.1073/pnas.0506580102.
#' @author Amanda Hiser
#' @examples
#' \dontrun{
#' # Prepare ranks
#' ranks <- t_tests$t_stat
#' names(ranks) <- t_tests$EntrezGeneSymbol
#' ranks <- sort(ranks, decreasing = TRUE)
#'
#' # Run GSEA
#' gsea_res <- somaPrGSEA(ranks = ranks)
#'
#' # Plot without leading edge annotation, using first row/pathway (by default)
#' plotES(gsea_results = gsea_res)
#'
#' # Plot with leading edge annotation for a specified pathway
#' plotES(x = 5, gsea_res, show_leading_edge = TRUE)
#' 
#' # Multiple plots can be generated with [lapply()]
#' top_paths <- gsea_res$results$pathway[1:3]
#' lapply(top_paths, function(x) plotES(gsea_results = gsea_res, 
#'                                      x = x, 
#'                                      show_leading_edge = TRUE))
#' }
#' @importFrom ggplot2 element_line element_text margin theme
#' @importFrom patchwork wrap_plots
#' @export
plotES <- function(x = 1,
                   gsea_results,
                   cust_path = NULL,
                   show_leading_edge = FALSE,
                   show_path_stats = TRUE,
                   enrichment_score_color = "green",
                   font_size = 12,
                   stats_color = "black",
                   leading_edge_color = "red",
                   leading_edge_size = 2,
                   leading_edge_alpha = 0.8) {

    stopifnot(
        "`gsea_results` must be the list output of somaPrGSEA()" = inherits(gsea_results, "list"),
        "`show_leading_edge` must be logical" = is.logical(show_leading_edge)
    )
    
    res_df <- gsea_results$results
    res_ranks <- sort(gsea_results$final_ranks, decreasing = TRUE) # To match what fgsea does internally
    
    # Retrieve pathway specified by x arg
    if ( is.numeric(x) ) {
        row_idx <- x # Get row index
    } else {
        row_idx <- which(res_df$pathway == x) # Get pathway name
        if ( length(row_idx) == 0 ) {
            stop("x '", x, "' not found in 'gsea_results'.")
        }
        row_idx <- row_idx[1]
    }
    
    # Extract leading edge genes & other relevant info
    le_feats <- res_df$leadingEdge[[row_idx]]
    pthwy_name <- res_df$pathway[[row_idx]]
    
    if ( !is.null(cust_path) ) {
        stopifnot(
            "`cust_path` must be a named list." = inherits(cust_path, "list") && !is.null(names(cust_path))
        )
        pathway <- unlist(cust_path[[pthwy_name]]) # cust_path must be named!
        pathway <- unique(pathway)
    } else {
        pid     <- res_df$pathway_id[[row_idx]]
        pathway <- pathway_map[pathway_map$pathway_id == pid, ]$gene_symbol
    }
    
    # Extract NES and ES values for annotation
    nes_val <- res_df$NES[[row_idx]]
    es_val <- res_df$ES[[row_idx]]
    n_pathway_genes <- length(pathway)
    n_le_genes <- ifelse(!is.null(le_feats), length(le_feats), 0L)
    
    if ( is.null(le_feats) ) {
        le_feats <- character(0)
        message("No leading edge genes found; none will be annotated on the plot.")
    }

    # Calculate running enrichment score data
    plot_df <- .calc_plot_data(res_ranks, pathway)

    # Mark leading edge genes if show_leading_edge is TRUE.
    # Need additional logic to handle multimers, otherwise they will be 
    # plotted multiple times
    if ( show_leading_edge && length(le_feats) > 0 ) {
        if ( es_val >= 0 ) {
            es_peak_idx <- which.max(plot_df$runningScore)
            is_candidate <- plot_df$gene %in% le_feats &
                            plot_df$position == 1 &
                            plot_df$x <= es_peak_idx
        } else {
            es_peak_idx <- which.min(plot_df$runningScore)
            is_candidate <- plot_df$gene %in% le_feats &
                            plot_df$position == 1 &
                            plot_df$x >= es_peak_idx
        }
        # If a gene name appears multiple times in the ranked list (e.g.
        # multimers), only mark first occurrence
        le_genes_seen <- plot_df$gene %in% le_feats & is_candidate
        is_dup <- duplicated(ifelse(le_genes_seen, plot_df$gene, NA)) & le_genes_seen
        plot_df$is_leading_edge <- is_candidate & !is_dup
    }

    p_res <- .create_es_plot(plot_df = plot_df,
                             color = enrichment_score_color,
                             font_size = font_size,
                             show_leading_edge  = show_leading_edge,
                             leading_edge_color = leading_edge_color,
                             leading_edge_size  = leading_edge_size,
                             leading_edge_alpha = leading_edge_alpha,
                             title = pthwy_name)

    p_hits <- .create_hits_plot(plot_df   = plot_df, 
                                font_size = font_size,
                                show_leading_edge  = show_leading_edge,
                                leading_edge_color = leading_edge_color)

    p_ranked <- .create_rank_plot(plot_df   = plot_df, 
                                  font_size = font_size,
                                  nes = nes_val,
                                  es = es_val,
                                  n_genes = n_pathway_genes,
                                  n_le = n_le_genes,
                                  show_path_stats = show_path_stats,
                                  stats_color = stats_color)

    # Combine selected subplots
    plotlist <- list(p_res, p_hits, p_ranked)
    n <- length(plotlist)

    # Add x-axis elements to the bottom plot
    plotlist[[n]] <- plotlist[[n]] +
        theme(axis.line.x  = element_line(),
              axis.ticks.x = element_line(),
              axis.text.x  = element_text())

    # Combine plots using patchwork
    combined <- patchwork::wrap_plots(plotlist, 
                                      ncol = 1, 
                                      heights = c(1.5, 0.3, 1))
    
    return(combined)
}


#' Calculate GSEA plot data (running enrichment scores)
#'
#' Calculates the running enrichment score for a gene set across a ranked
#' gene list, following the GSEA algorithm.
#'
#' @param ranks Named numeric vector of ranked genes (sorted by rank).
#' @param geneSet Character vector of genes in the gene set.
#' @param exponent Numeric. Weight exponent for the running sum calculation.
#'   Default is 1.
#' @return A data.frame with columns:
#'   \item{x}{Position in ranked list (1 to n)}
#'   \item{runningScore}{Cumulative enrichment score at each position}
#'   \item{position}{1 if gene is in the gene set, 0 otherwise}
#'   \item{gene}{Gene identifier}
#'   \item{geneList}{The ranking statistic value}
#'
#' @noRd
.calc_plot_data <- function(ranks, geneSet, exponent = 1) {
    
    n <- length(ranks) # Total number of features
    gene_names <- names(ranks)

    # Identify which genes are in the gene set
    hits <- gene_names %in% geneSet

    # Calculate the weighted scores for hits
    # Use absolute value of statistic raised to the exponent
    hit_scores <- abs(ranks)^exponent * hits

    # Normalization factors
    nr <- sum(hit_scores)
    n_miss <- n - sum(hits)

    # Calculate running score
    if (nr > 0 && n_miss > 0) {
        hit_contribution <- hit_scores / nr
        miss_contribution <- (!hits) / n_miss
        score_changes <- hit_contribution - miss_contribution
    } else if (nr == 0) {
        score_changes <- rep(0, n) # No hits in the gene set
    } else {
        hit_contribution <- hit_scores / nr
        score_changes <- hit_contribution # All genes are hits
    }

    running_score <- cumsum(score_changes)

    # Create output data frame
    data.frame(x = seq_len(n),
               runningScore = running_score,
               position = as.integer(hits),
               gene = gene_names,
               geneList = as.numeric(ranks))
}


#' Create the running enrichment score subplot
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point annotate
#' @importFrom ggplot2 labs theme element_text element_blank element_line margin scale_x_continuous
#' @noRd
.create_es_plot <- function(plot_df, color, font_size,
                            show_leading_edge, leading_edge_color, leading_edge_size,
                            leading_edge_alpha, title) {

    p <- ggplot(plot_df, aes(x = x)) +
        geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
        geom_line(aes(y = runningScore), color = color, linewidth = 1) +
        labs(x = NULL, y = "Running Enrichment Score", title = title) +
        theme(panel.background   = element_blank(),
              plot.background    = element_blank(),
              panel.border       = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line.y        = element_line(colour = "black"),
              axis.line.x        = element_blank(),
              axis.text.x        = element_blank(),
              axis.ticks.x       = element_blank(),
              plot.margin        = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm"),
              plot.title         = element_text(hjust = 0.5)) +
        scale_x_continuous(expand = c(0, 0))

    # Add leading edge gene points
    if ( show_leading_edge ) {
        le_data <- plot_df[plot_df$is_leading_edge, ]

        if ( nrow(le_data) > 0 ) {
            p <- p + geom_point(data = le_data,
                                aes(x = x, y = runningScore),
                                color = leading_edge_color,
                                size  = leading_edge_size,
                                alpha = leading_edge_alpha,
                                shape = 19)

                # Position legend in upper right of plot area
                max_x <- max(plot_df$x)
                score_range <- range(plot_df$runningScore)
                legend_y <- score_range[2] - diff(score_range) * 0.1

                # Add 'Leading Edge' legend w/ point symbol
                p <- p + 
                    annotate("point",
                             x     = max_x * 0.82,
                             y     = legend_y,
                             color = leading_edge_color,
                             size  = leading_edge_size,
                             shape = 19) +
                    annotate("text",
                             x     = max_x * 0.84,
                             y     = legend_y,
                             label = "Leading Edge",
                             hjust = 0,
                             size  = font_size / 3)
        }
    }

    return(p)
}


#' Create the gene hit positions subplot
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_linerange labs
#' @importFrom ggplot2 theme element_blank element_line margin
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_color_identity
#' @noRd
.create_hits_plot <- function(plot_df, font_size,
                              show_leading_edge = FALSE,
                              leading_edge_color = "red") {

    hits_data <- plot_df[plot_df$position == 1, ] # Filter to only hits

    # Assign a color column so LE hits match the ES curve annotation
    if ( show_leading_edge && "is_leading_edge" %in% names(hits_data) ) {
        hits_data$hit_color <- ifelse(hits_data$is_leading_edge,
                                      leading_edge_color, "black")
    } else {
        hits_data$hit_color <- "black"
    }

    ggplot(hits_data, aes(x = x)) +
        geom_linerange(aes(ymin = 0, ymax = 1, color = hit_color),
                       linewidth = 0.3) +
        scale_color_identity() +
        labs(x = NULL, y = NULL) +
        theme(legend.position   = "none",
              panel.background   = element_blank(),
              plot.background    = element_blank(),
              panel.border       = element_blank(),
              panel.grid.major   = element_blank(),
              panel.grid.minor   = element_blank(),
              axis.line.y        = element_line(colour = "black"),
              axis.line.x        = element_blank(),
              plot.margin        = margin(t = -.1, b = 0, unit = "cm"),
              axis.ticks         = element_blank(),
              axis.text          = element_blank()) +
        scale_x_continuous(limits = c(1, nrow(plot_df)), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
}


#' Create the ranked list metric subplot
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_segment geom_hline annotate
#' @importFrom ggplot2 labs theme element_blank element_line margin scale_x_continuous
#' @noRd
.create_rank_plot <- function(plot_df, font_size, stats_color,
                              nes, es, n_genes, n_le,
                              show_path_stats = TRUE) {

    # Filter to only hits (for plotting segments only)
    hits_data <- plot_df[plot_df$position == 1, ]

    p <- ggplot(plot_df, aes(x = x)) +
        geom_segment(data = hits_data, aes(x = x, xend = x, 
                                           y = geneList, yend = 0),
                     color = "grey50") +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
        theme(panel.background   = element_blank(),
              plot.background    = element_blank(),
              panel.border       = element_blank(),
              panel.grid.major   = element_blank(),
              panel.grid.minor   = element_blank(),
              axis.line.x        = element_line(colour = "black"),
              axis.line.y        = element_line(colour = "black"),
              plot.margin        = margin(t = -.1, r = .2, b = .2, l = .2, 
                                         unit = "cm")) +
        scale_x_continuous(limits = c(1, nrow(plot_df)), expand = c(0, 0))
    
    # Add annotations for NES, ES, pathway size, and leading edge count
    if ( show_path_stats ) {
        max_x <- max(plot_df$x)
        metric_range <- range(plot_df$geneList)
    
        # Create annotation label
        annot_label <- paste0("ES = ", round(es, 3), "\n",
                              "NES = ", round(nes, 3), "\n",
                              "Pathway genes (N) = ", n_genes, "\n",
                              "Leading edge (N) = ", n_le)
    
        # Position in bottom right corner
        annot_y <- metric_range[1] + diff(metric_range) * 0.05
    
        p <- p + annotate("text",
                          x = max_x * 0.98,
                          y = annot_y,
                          label = annot_label,
                          hjust = 1,
                          vjust = 0,
                          size = font_size / 3,
                          color = stats_color)
    }
    
    return(p)
}
