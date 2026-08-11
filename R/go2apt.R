#' Retrieve SomaScan Analytes in a Given GO Term
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
    
    if ( is.null(col_meta_df) ) {
        if ( verbose ) {
            message("`col_meta_df` not provided, using 11K annotations by default.")
        }
        col_meta_df <- SomaDataIO::getAnalyteInfo(example_data_11k)
    }
    
    go_map <- pathway_map[pathway_map$group_code %in% c("bp", "mf"), ]
    
    if ( !x %in% go_map$pathway_id ) {
        err_msg <- paste0("The provided GO term '", x, "'", " was not found.")
        stop(err_msg)
    }
    
    genes <- go_map[go_map$pathway_id == x, ]$gene_symbol
    apts <- gene2apt(genes, col_meta_df, collapse = FALSE)
    apts <- unique(apts)
    
    return(apts[!is.na(apts)])
}
