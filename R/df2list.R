#' Convert a Data Frame to a Named List of Vectors
#'
#' @param df Data frame to be converted.
#' @param name_col Character. Name of the column in `df` to be used for list names.
#' @param value_col Character. Name of the column in `df` to be used as list
#'   elements (vectors).
#' @param clean_names Logical. Should list names be modified to remove
#'   whitespaces and non-aphanumeric strings (excluding underscores)? See
#'   [SomaDataIO::cleanNames()] for more details.
#' @returns A list of vectors.
#' @author Amanda Hiser
#' @examples
#' df <- SomaDataIO::rn2col(SomaDataIO::example_data, "Rowname")
#' df2list(df, name_col = "SampleType", value_col = "Rowname")
#'
#' # Use clean_names argument to make R-compatible list names
#' df <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
#' head(df$TargetFullName)
#' df2list(df, name_col = "TargetFullName", value_col = "AptName",
#'         clean_names = TRUE) |>
#' head()
#' @importFrom SomaDataIO cleanNames
#' @export
df2list <- function(df, name_col, value_col, clean_names = FALSE) {
  stopifnot(inherits(df, "data.frame"))
  stopifnot(
    "`name_col` must be a column name in `df`" = name_col %in% colnames(df),
    "`value_col` must be a column name in `df`" = value_col %in% colnames(df)
  )
  if ( clean_names ) {
    df[[name_col]] <- SomaDataIO::cleanNames(df[[name_col]])
  }
  df <- df[, c(name_col, value_col)]

  # Set as factor to retain original order in final list
  name_var <- factor(df[[name_col]], levels = unique(df[[name_col]]))

  split(df[[value_col]], name_var)
}
