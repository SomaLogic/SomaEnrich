#' Convert Gene Set into SomaScan ID-Based Feature Set
#'
#' @description
#' Converts gene pathways to SomaScan identifiers. All analytes that map to
#' a given gene are retained; if multiple analytes map to the same gene, all
#' will be present in the output.
#'
#' @param path_list Named list of character vectors containing genes.
#' @param col_meta_df Data frame of column metadata, likely created by
#'   [SomaDataIO::getAnalyteInfo()]. This object should contain SomaScan menu 
#'   information that can be used to map between `AptNames` and Entrez gene 
#'   identifiers. At minimum, it must contain columns named "AptName" 
#'   (containing SomaScan identifiers in `AptName` format) and either 
#'   "EntrezGeneSymbol" or "EntrezGeneID". 
#' @param id_type Character. Type of gene identifier used in `path_list`. 
#'   Options are "EntrezGeneSymbol" or "EntrezGeneID". Default is 
#'   "EntrezGeneSymbol".
#' @param verbose Logical. Show status messages and (if applicable) a progress 
#'   bar?
#' @returns List of named character vectors (names = genes, values = AptNames).
#'   When multiple analytes map to the same gene, all are retained as separate
#'   elements with the same name.
#' @author Amanda Hiser
#' @examples
#' genePath2aptPath(path_list = list(Example = ex_gene_pathway), 
#'                 col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k))
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
genePath2aptPath <- function(path_list, 
                             col_meta_df,
                             id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                             verbose = interactive()) {
    
    id_type <- match.arg(id_type)
    id_col  <- .match_col(col_meta_df, id_type)
    
    # Create mapping for genes -> AptNames.
    # Build once, outside of lapply call
    id_map <- .make_id_map(col_meta_df, 
                           gene_col = id_col, 
                           apt_col  = "AptName")
    
    n_pth <- length(path_list) # Use for progress messaging
    
    if (verbose && n_pth > 1) {
        message("Converting ", n_pth, " gene-based pathways to AptName format...")
        pb <- txtProgressBar(min   = 0, 
                             max   = n_pth, 
                             style = 3)
    }
    
    # Convert each path from genes -> AptNames
    conv_paths <- lapply(seq_along(path_list), function(i) {
        if (verbose && n_pth > 1) setTxtProgressBar(pb, i)
        .convGenePath(pathway      = path_list[[i]], 
                      gene_apt_map = id_map, 
                      id_type      = id_col)
    })
    
    if (verbose && n_pth > 1) {
        close(pb)
        message("Done!")
    }
    
    names(conv_paths) <- names(path_list)
    conv_paths
}


.convGenePath <- function(pathway, gene_apt_map, id_type) {
    
    genes <- unique(unname(unlist(pathway)))
    
    # Build pathway-specific lookup table
    lookup_df <- gene_apt_map[gene_apt_map[[id_type]] %in% genes, ]
    lookup_df <- lookup_df[!is.na(lookup_df$AptName), ]
    
    if (nrow(lookup_df) == 0) {
        return(setNames(character(0), character(0)))
    }
    
    final_path <- lookup_df$AptName
    names(final_path) <- lookup_df[[id_type]]
    
    return(final_path)
}

