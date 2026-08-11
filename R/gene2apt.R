#' Convert Gene Identifiers to SomaScan Analytes (v2)
#' 
#' Takes a set of gene identifiers and returns associated SomaScan
#' analytes, given a `col_meta` object returned from
#' [SomaDataIO::getAnalyteInfo()].
#'
#' @inheritParams apt2gene 
#' @param x Character. A vector of gene identifiers.
#' @param id_type Character. Type of gene identifier provided in `x`.
#'   Options are "EntrezGeneSymbol" or "EntrezGeneID". Default is
#'   "EntrezGeneSymbol".
#' @param sep Character used to separate collapsed identifiers. Default is "|".
#' @returns A character vector of SomaScan identifiers (in `AptName` format).
#' @author Amanda Hiser, Alex Poole
#' @examples
#' anno <- SomaDataIO::getAnalyteInfo(example_data_11k)
#' genes <- withr::with_seed(123, sample(anno$EntrezGeneSymbol, 5))
#' gene2apt(x = genes, col_meta_df = anno)
#' 
#' # By default, one string for each gene will be returned.
#' # Genes targeted by >1 analytes will have analytes collapsed via `sep`
#' gene2apt(x = "IL12B", col_meta_df = anno)
#' 
#' # Use `collapse = FALSE` to return all values individually
#' gene2apt(x = "IL12B", col_meta_df = anno, collapse = FALSE)
#' @importFrom SomaDataIO getAnalyteInfo
#' @importFrom tidyr separate_rows
#' @importFrom stats na.omit
#' @export
gene2apt <- function(x,
                     col_meta_df = NULL,
                     id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                     collapse = TRUE,
                     sep = "|") {
    
    id_type <- match.arg(id_type)
    
    if ( inherits(x, "list") ) {
        stop("`x` must be a vector, not a list", call. = FALSE)
    }
    
    if ( is.null(col_meta_df) ) {
        col_meta_df <- SomaDataIO::getAnalyteInfo(example_data_11k)
    } else {
        stopifnot("`col_meta_df` must be a data frame or tibble." = inherits(col_meta_df, "data.frame"))
    }
    
    id_col <- .match_col(col_meta_df, id_type)

    # Clean input by converting blanks to NA
    x[x == ""] <- NA_character_
    
    id_map <- .make_id_map(x = x,
                           col_meta_df = col_meta_df,
                           gene_col = id_col,
                           apt_col = "AptName")
    
    if ( collapse ) {
        collapsed <- .collapse_by_gene(id_map,
                                       apt_col = "AptName",
                                       gene_col = id_col,
                                       sep = sep)
        # Align results back to input (preserving order, duplicates, and NA for unmatched)
        x_clean <- x[!is.na(x)]
        unname(collapsed[match(x_clean, names(collapsed))])
    } else {
        unique(stats::na.omit(id_map[["AptName"]]))
    }
}
