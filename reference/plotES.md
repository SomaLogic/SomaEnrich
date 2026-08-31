# GSEA Enrichment Plot with Leading Edge Annotation

Creates a GSEA enrichment plot, with the option to annotate leading edge
features on the running enrichment score curve.

## Usage

``` r
plotES(
  x = 1,
  gsea_results,
  cust_path = NULL,
  show_leading_edge = FALSE,
  show_path_stats = TRUE,
  enrichment_score_color = "green",
  font_size = 12,
  stats_color = "black",
  leading_edge_color = "red",
  leading_edge_size = 2,
  leading_edge_alpha = 0.8
)
```

## Arguments

- x:

  Character or integer. The pathway name *or* row index of
  `gsea_results` to use for creating the plot and extracting leading
  edge features. Default is 1 (first row). If a character, must match a
  pathway name in the "pathway_name" column of `gsea_results`.

- gsea_results:

  A `data.frame` of GSEA results from
  [`somaPrGSEA()`](https://somalogic.github.io/SomaEnrich/reference/somaPrGSEA.md).
  The `leadingEdge` column will be used to identify and annotate leading
  edge features on the plot. The `resource_id` column will be used to
  retrieve pathway members. If a custom pathway was used for GSEA, the
  pathway must be provided to `cust_path`.

- cust_path:

  Optional. A vector of gene identifiers (gene symbols or Entrez gene
  IDs) or SomaScan `AptNames` representing custom gene sets or
  functional groups of interest. Should be provided if the `cust_paths`
  argument was used to generate results via
  [`somaPrGSEA()`](https://somalogic.github.io/SomaEnrich/reference/somaPrGSEA.md).

- show_leading_edge:

  Logical. Should leading edge features be highlighted on the plot? If
  TRUE, leading edge features are extracted from `gsea_results` and
  annotated with colored dots. Default is TRUE.

- show_path_stats:

  Logical. Should GSEA result statistics (ES, NES, and path size) be
  annotated on the plot? Default is TRUE.

- enrichment_score_color:

  Character. Color of the running enrichment score line. Default is
  "green".

- font_size:

  Numeric. Base font size for the plot. Default is 12.

- stats_color:

  Character. Color of the stats annotation from `show_path_stats`.

- leading_edge_color:

  Character. Color of the leading edge gene dots. Default is "red".

- leading_edge_size:

  Numeric. Size of the leading edge gene dots. Default is 2.

- leading_edge_alpha:

  Numeric. Transparency of the leading edge dots. Default is 0.8.

## Value

A ggplot object (combined via
[patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html))
displaying the GSEA enrichment plot with optional leading edge gene
annotation.

## Details

This function creates a plot with three panels:

1.  **Enrichment Score Panel**: Shows the running enrichment score with
    leading edge genes optionally highlighted as colored dots.

2.  **Gene Hits Panel**: Shows vertical lines indicating where genes in
    the gene set appear in the ranked list.

3.  **Ranked List Metric Panel**: Shows the distribution of the ranking
    metric across all genes.

The leading edge genes are the core subset that contribute most to the
enrichment signal, appearing at or before the point where the running
enrichment score reaches its maximum deviation from zero. For a positive
ES, the leading edge subset is the set of members that appear in the
ranked list *prior* to the peak score. For a negative ES, the leading
edge subset will appear subsequent to the peak score.

## References

Subramanian, Tamayo, et al. Gene set enrichment analysis: A
knowledge-based approach for interpreting genome-wide expression
profiles. Proc Natl Acad Sci USA. 102(43):15545-50 (2005).
https://doi.org/10.1073/pnas.0506580102.

## Author

Amanda Hiser

## Examples

``` r
if (FALSE) { # \dontrun{
# Prepare ranks
ranks <- t_tests$t_stat
names(ranks) <- t_tests$EntrezGeneSymbol
ranks <- sort(ranks, decreasing = TRUE)

# Run GSEA
gsea_res <- somaPrGSEA(ranks = ranks)

# Plot without leading edge annotation, using first row/pathway (by default)
plotES(gsea_results = gsea_res)

# Plot with leading edge annotation for a specified pathway
plotES(x = 5, gsea_res, show_leading_edge = TRUE)

# Multiple plots can be generated with [lapply()]
top_paths <- gsea_res$results$pathway[1:3]
lapply(top_paths, function(x) plotES(gsea_results = gsea_res, 
                                     x = x, 
                                     show_leading_edge = TRUE))
} # }
```
