#' Retrieve SomaScan Analytes in a Given GO Term
#' 
#' @details
#' A wrapper around [path2apt()] specifically for GO terms.
#' 
#' @inheritParams apt2gene
#' @param x Character. A single GO term ID. Must contain "GO:" prefix.
#' @returns Character vector of SomaScan analyte identifiers in `AptName`
#'   format. All analytes associated with each gene are returned. For more
#'   information about SomaScan identifiers and their formats, please see
#'   [SomaDataIO::SeqId].
#' @importFrom SomaDataIO getAnalyteInfo
#' @export
#' @examples
#' go2apt("GO:0019319") # Returns vector of AptNames
#' 
#' # Apply over a vector of GO terms
#' sapply(c("GO:2001256", "GO:0000462", "GO:2001259"), go2apt, USE.NAMES = TRUE)
go2apt <- function(x,
                   col_meta_df = NULL,
                   verbose = interactive()) {
    
    if ( length(x) > 1 ) {
        stop("`go2apt()` accepts only 1 GO term as input.", call. = FALSE)
    }
    
    go_ids <- pathway_map$pathway_id[pathway_map$group_code %in% c("bp", "mf")]
    
    if ( !x %in% go_ids ) {
        err_msg <- paste0("The provided GO term '", x, "'", " was not found.")
        stop(err_msg, call. = FALSE)
    }
    
    # Ensures a vector is returned, instead of a list
    results <- path2apt(x, col_meta_df = col_meta_df, verbose = verbose)[[1L]]
    
    return(results)
}
