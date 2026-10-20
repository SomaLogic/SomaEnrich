#' Create Nominal Significance Table
#'
#' Builds a significance table across combinations of nominal p-value and
#' fold-change cutoffs, typically used when FDR-adjusted p-values do not
#' reach the desired alpha. This function is not meant to be a substitute for
#' formal multiple-testing correction. Rather, it is useful for small or pilot
#' studies, power planning, or assessing how sensitive results are to
#' chosen cutoffs.
#'
#' @details
#' This function counts how many analytes in a contrast results table
#' pass each combination of p-value and fold-change cutoffs. For every
#' cutoff pair, it returns the number of analytes that meet both criteria
#' in each contrast. This function is best used to explore how many
#' analytes would be called significant under different thresholds,
#' especially when FDR-adjusted results are sparse.
#'
#' \describe{
#'   \item{When to use:}{After running pairwise contrasts, when few or no
#'   analytes pass FDR, and you want to see how many would pass at common
#'   nominal p-value and effect size thresholds (e.g. p ≤ 0.05 and
#'   |FC| ≥ 0.5).}
#'   \item{When _not_ to use:}{Do not treat these counts as FDR-controlled
#'   "significant" findings for publication or confirmatory claims. Prefer
#'   FDR-adjusted results when sample size and power support them.}
#' }
#'
#' Rows with a missing (`NA`) fold change are excluded before counting.
#' When `group_col` is used, missing group values are counted as their own
#' `NA` group, and every group appears at every cutoff combination
#' (with a `Count` of 0 when no analytes in that group pass).
#'
#' @param results A `data.frame` of pairwise contrast results, like those
#'   generated from the R package `emmeans`. Must contain the columns
#'   used for the `stat_col`, `fc_col`, and `group_col` arguments.
#' @param stat_col Character. Column name for the p-value or FDR-adjusted
#'   p-value used for significance filtering. Default is "p.value".
#' @param fc_col Character. Column name for the fold change (or effect size)
#'   used for magnitude filtering. Default is "estimate" (from `emmeans`
#'   output).
#' @param group_col Character or `NULL`. Quoted column name identifying the
#'   contrast (grouping variable). If a column name is provided, counts are
#'   computed separately within each level of that column. The column's
#'   type is preserved, and groups are ordered by factor level (or sorted,
#'   for other types). If `NULL`, counts are computed across all rows
#'   without grouping. Default is `NULL`.
#' @param sig.levels Numeric vector of p-value cutoffs to evaluate.
#'   Default is `c(0.1, 0.05, 0.01)`.
#' @param fc.levels Numeric vector of fold change cutoffs to evaluate.
#'   Default is `c(0.25, 0.5, 0.75)`.
#' @returns A `data.frame` with one row for each combination of `sig.levels`
#'   and `fc.levels`. When `group_col` is used, each combination has one row
#'   per group. The columns are:
#'   \item{group_col}{only present if `group_col` is not `NULL`}
#'   \item{cutoff}{the p-value cutoff}
#'   \item{fc}{the fold change cutoff}
#'   \item{Count}{the number of analytes with a p-value at or below
#'     `cutoff` and an absolute fold change at or above `fc`}
#' @examples
#' # Synthesize emmeans-style pairwise contrast results
#' contrasts <- withr::with_seed(123, {
#'   data.frame(AptName = rep(sprintf("seq.%d.%02d", sample(10000:99999, 6),
#'                           sample(0:99, 6)), each = 2),
#'             time_contrast = rep(c("Week4 - Week0", "Week12 - Week0"), 6),
#'             estimate = c(0.10, 0.55, 0.80, 0.30, -0.60, 0.20,
#'                          0.45, -0.90, 0.05, 0.70, -0.35, 1.10),
#'             p.value  = c(0.40, 0.02, 0.001, 0.08, 0.03, 0.50,
#'                          0.04, 0.0005, 0.70, 0.01, 0.12, 0.002))
#' })
#'
#' # No grouping
#' makeNominalTable(contrasts)
#'
#' # Group by time contrast
#' makeNominalTable(contrasts, group_col = "time_contrast")
#'
#' # Restrict to a single p-value/fold-change combination
#' makeNominalTable(contrasts,
#'                  sig.levels = 0.05,
#'                  fc.levels  = 0.5)
#' @export
makeNominalTable <- function(results,
                             stat_col   = "p.value",
                             fc_col     = "estimate",
                             group_col  = NULL,
                             sig.levels = c(0.1, 0.05, 0.01),
                             fc.levels  = c(0.25, 0.5, 0.75)) {

  # Check to make sure input has necessary cols
  missing_cols <- setdiff(c(stat_col, fc_col, group_col), names(results))
  
  if ( length(missing_cols) > 0L ) {
    stop("`results` is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if ( !is.null(group_col) ) {
    stopifnot("`group_col` must be a single character string." =
                is.character(group_col) && length(group_col) == 1L)

    # The output already has columns with these names, so a grouping column
    # with the same name would overwrite them
    reserved <- c("cutoff", "fc", "Count")
    
    if ( group_col %in% reserved ) {
      stop("`group_col` cannot be one of the reserved output column names: ",
           paste(reserved, collapse = ", "),
           ". Please rename that column in `results` first.",
           call. = FALSE)
    }
  }

  stopifnot("`sig.levels` must be a numeric vector." = is.numeric(sig.levels),
            "`fc.levels` must be a numeric vector."  = is.numeric(fc.levels))

  # Pull out the p-values and absolute fold changes used for filtering.
  # Rows missing either value are never counted.
  p      <- results[[stat_col]]
  abs_fc <- abs(results[[fc_col]])
  valid  <- !is.na(p) & !is.na(abs_fc)

  # Without grouping, treat every row as part of one group, so the same
  # counting code below works whether or not `group_col` is given.
  if ( is.null(group_col) ) {
    groups    <- NULL
    group_idx <- rep(1L, nrow(results))
  } else {
    # With grouping, list each group once (factors keep their level order,
    # other types are sorted, and NA goes last)
    groups    <- sort(unique(results[[group_col]]), 
                      method = "radix", 
                      na.last = TRUE)
    # Records which group each row belongs to
    group_idx <- match(results[[group_col]], groups)
  }

  # One row per cutoff pair and group, sorted by cutoff, then fc, then group,
  # so groups with no passing analytes still get a Count of 0
  out <- expand.grid(group  = seq_len(max(1L, length(groups))),
                     fc     = fc.levels,
                     cutoff = sig.levels,
                     KEEP.OUT.ATTRS = FALSE)

  # For each row of the table, count the analytes in that group that pass
  # both the p-value cutoff and the fold-change cutoff
  out$Count <- vapply(seq_len(nrow(out)), function(i) {
    sum(valid & group_idx == out$group[i] &
          p <= out$cutoff[i] & abs_fc >= out$fc[i])
  }, integer(1))

  # Replace group numbers with actual group labels
  if ( !is.null(group_col) ) {
      out[[group_col]] <- groups[out$group]
  }

  out[c(group_col, "cutoff", "fc", "Count")]
}
