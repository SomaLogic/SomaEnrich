#' Check For Pre-Sorting
#' 
#' Checks to see if a vector has been sorted in either ascending or descending
#' order.
#' @noRd
.is.sorted <- function(x, ...) {
  !is.unsorted(x, na.rm = TRUE, ...) | !is.unsorted(rev(x), na.rm = TRUE, ...)
}


#' Split Delimited IDs
#' 
#' Used to split up gene IDs delimited with "|" (or other characters) from the 
#' SomaScan menu. Used for pre-processing in both GSEA and ORA.
#' @noRd
.splitIDs <- function(ids, sep = "\\,|\\s+|\\|") {
    
    if ( !is.null(names(ids)) ) {
        # Splits names and duplicates values for the new split names
        split_names <- strsplit(names(ids), sep)
        ids_split <- setNames(rep(ids, lengths(split_names)), unlist(split_names))
    } else {
        ids_split <- unlist(strsplit(ids, sep))
    }
    
    return(ids_split)
}


# Borrow from cli for internal use, without explicitly importing
# the package in the NAMESPACE file
.todo <- cli::cli_alert
.done <- cli::cli_alert_success
.oops <- cli::cli_alert_danger
.warn <- cli::cli_warn
.info <- cli::cli_alert_info
.rule <- cli::rule


#' Match IPP and Array-Derived ADAT Columns
#'
#' ADATs produced via Illumina Protein Prep (IPP) may have different field
#' names than ADATs produced by SomaLogic's array-based readout. This function
#' allows more flexibility when matching column metadata names, ensuring that
#' fields from both IPP NGS-readout ADATs _and_ array-readout ADATs can be
#' matched.
#' 
#' Returns the actual column name in `df` that matches `col_name`. First, tries
#' an exact match; if that fails, strips non-alphanumeric 
#' characters and compares lowercased strings to find a fuzzy match.
#' @param df A data frame.
#' @param col_name Character. The column name to match (e.g. `"EntrezGeneSymbol"`).
#' @return The matched column name as it appears in `names(df)`.
#' @keywords internal
#' @noRd
.match_col <- function(df, col_name, arg = "id_type") {
    # Return early in cases of exact matches
    if ( col_name %in% names(df) ) {
      return(col_name)
    }
    # Normalize col names if exact match isn't found
    norm <- function(x) gsub("[^[:alnum:]]", "", tolower(x))
    hit  <- names(df)[norm(names(df)) == norm(col_name)]
    
    if ( length(hit) == 1L ) {
      return(hit)
    }
    
    if ( length(hit)  > 1L ) {
        stop("Ambiguous column match for '", col_name, "': ",
             paste(hit, collapse = ", "), call. = FALSE)
    }
    
    stop("No column matching '", col_name, "' found in the provided data frame. ",
         "Check that `", arg, "` corresponds to an existing column.", call. = FALSE)
}


###### Identifier conversion sub-functions ##########

#' Build Gene to Aptamer Mapping Table
#' 
#' @param col_meta_df Data frame of column metadata from [SomaDataIO::getAnalyteInfo()].
#' @param x Optional. List of genes or analytes of interest. Will be used to filter 
#'   `col_meta_df`. If not provided, all features in `col_meta_df` are considered.
#' @param gene_col Character. Name of the column containing gene identifiers.
#' @param apt_col Character. Name of the column containing AptNames.
#' @return A data frame with one row per gene-aptamer pair. Empty/NA genes are
#'   removed from the reference map.
#' @keywords internal
#' @importFrom SomaDataIO getAnalyteInfo
#' @importFrom tidyr separate_rows
#' @noRd
.make_id_map <- function(col_meta_df, 
                         x = NULL,
                         gene_col = "EntrezGeneSymbol",
                         apt_col = "AptName") {
    
    # Determine what identifier is being provided
    id_col <- if ( is.null(x) ) { 
        apt_col # Use AptNames by default
    } else {
        ifelse(all(SomaDataIO::is.AptName(x)), apt_col, gene_col)
    }
    
    # Expand multi-gene entries into separate rows.
    # TODO: Find an elegant way to achieve this w/o tidyr
    exp_df <- col_meta_df[, c(apt_col, gene_col)] |>
        tidyr::separate_rows(!!gene_col, sep = "\\||\\,|\\s+") |>
        as.data.frame()
    
    if ( !is.null(x) ) {
        exp_df <- exp_df[exp_df[[id_col]] %in% x, ]
    }
    
    # Remove empty and NA gene entries
    keep <- exp_df[[gene_col]] != "" & !is.na(exp_df[[gene_col]])
    id_map <- exp_df[keep, ]
    id_map <- unique(id_map) # Remove duplicate rows, if any
    
    return(id_map)
}


#' Collapse Analytes by Input Gene
#' 
#' @param id_map Data frame from `.make_id_map()`.
#' @param apt_col Character. Name of the column containing analytes.
#' @param gene_col Character. Name of the column containing genes.
#' @param sep Character. Separator for collapsing.
#' @return Character vector, one element per unique input gene. Returns NA
#'   for genes with no matches.
#' @keywords internal
#' @importFrom stats na.omit
#' @noRd
.collapse_by_gene <- function(id_map, 
                              apt_col = "AptName", 
                              gene_col = "EntrezGeneSymbol",
                              sep = "|") {
    
    # Split by gene ID, each list name is a gene & elements include AptNames
    split_df <- split(id_map, id_map[[gene_col]])
    
    collapsed_vec <- vapply(split_df, function(df) {
        apts <- unique(stats::na.omit(df[[apt_col]]))
        
        if ( length(apts) == 0 ) {
            NA_character_ # Return NA if no AptNames found for the gene
        } else {
            paste(apts, collapse = sep) # If multiple found, collapse into 1 string 
        }
    }, character(1), USE.NAMES = TRUE)
    
    return(collapsed_vec)
}


#' Select Best Value Based on Metric
#' 
#' When multiple values map to the same protein/gene identifier, retains a 
#' single value to represent the gene based on a pre-defined selection method.
#' 
#' @param ranks A named numeric vector.
#' @param resolve_method Character. Selection method for choosing a single 
#'   representative value for a gene when non-1:1 mapping exists between 
#'   analytes and genes. Choose one of:
#'   \describe{
#'     \item{"abs"}{(Default) Selects the entry with the largest absolute value per
#'       gene. Recommended for signed statistics (e.g. t-statistic, log fold change)
#'       in bidirectional analyses — unlike \code{"rank"} or \code{"max"}, it 
#'       preserves strong signal symmetrically at both the top and bottom of 
#'       the ranking vector.}
#'     \item{"min"}{Selects the entry with the minimum numeric value, including
#'       negative values. Use for metrics where smaller values are better 
#'       (e.g., p-values) or when smaller/negative signed values are important.}
#'     \item{"max"}{Selects the entry with the maximum numeric value. Use for
#'       metrics where larger values are preferred (e.g. positive one-sided tests).}
#'     \item{"rank"}{Selects the entry with the smallest positional
#'       index (i.e. first occurrence in the input vector). Use when input is
#'       pre-sorted by one or more variables, and the "best" entry is the one closest 
#'       to position 1. This means that the selected entry depends on the
#'       order of the input, so `ranks` must be **pre-sorted** to place values
#'       of highest importance first.}
#'   }
#' @return List containing:
#'    \item{ranks}{Ranking vector used to select features. May differ from
#'                 input ranks if heterodimers are present, as they will be
#'                 split into individual components.}
#'    \item{id_list}{List of features. Elements are values from the `ranks`
#'                   vector.}
#'    \item{keep_list}{List of features, with only one value per feature. 
#'                     List elements are the index of the retained value (based
#'                     on "method") from the original `ranks` vector.}
#' @keywords internal
#' @importFrom rlang warn
.select_best <- function(ranks,
                         resolve_method = c("abs", "min", "max", "rank"),
                         call = rlang::caller_env()) {
    
    resolve_method <- match.arg(resolve_method)
    
    if ( resolve_method == "rank" && !.is.sorted(ranks) ) {
      rlang::warn(c("The `ranks` vector is not monotonically sorted.",
                  "When `resolve_method = 'rank'`, the first occurrence of each entry is retained, so the input order determines which value is kept.",
                  "If this is intentional (e.g. the data is sorted by a composite metric), this warning can be ignored."), 
                  call = call)
    }
    
    rank_df <- data.frame(id    = names(ranks),
                          value = unname(ranks),
                          order = seq_along(ranks))
    
    # Split data frame by ID to create list of data frames, grouped by ID
    df_list <- split(rank_df, rank_df$id)
    
    # Create list of vectors, with indices (order) as list element(s)
    id_list <- lapply(df_list, `[[`, "order")
    
    # Select "best" index for each group based on method. Duplicates that don't
    # meet selection criteria are dropped
    keep_list <- lapply(df_list, function(df) {
        idx <- switch(resolve_method,
                      rank = which.min(df$order),
                      min  = which.min(df$value),
                      max  = which.max(df$value),
                      abs  = which.max(abs(df$value))
        )
        df$order[idx] # Always returns the INDEX (order), regardless of method
    })
    
    return(list(ranks     = ranks,     # The cleaned ranks used in this funct
                id_list   = id_list,   # ALL indices, grouped by ID
                keep_list = keep_list)) # Index to retain for each ID
}
