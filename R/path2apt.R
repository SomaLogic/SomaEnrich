#' Retrieve SomaScan Analytes in Pathways from `pathway_map`
#'
#' Takes one or more pathway identifiers from `pathway_map$pathway_id` (any
#' resource, e.g. GO, MSigDB Hallmark, etc.) and returns the SomaScan analytes
#' that target the genes in each pathway. All analytes associated with each
#' gene are returned.
#'
#' @inheritParams apt2gene
#' @param x Character. A vector of pathway identifiers found in
#'   `pathway_map$pathway_id` (e.g. "GO:0019319", "M5890").
#' @param id_type Character. Gene identifier used to match pathway genes to
#'   analytes in `col_meta_df`; this does not change the output, which is
#'   always `AptNames`. Options are "EntrezGeneSymbol" or "EntrezGeneID".
#'   Entrez Gene IDs are robust to gene symbol changes. Default is
#'   "EntrezGeneSymbol".
#' @returns A named list of character vectors, one element per unique pathway
#'   in `x` (in the order provided). Each vector contains unique SomaScan
#'   analyte identifiers in `AptName` format. Pathways with no matching
#'   analytes in `col_meta_df` return `character(0)`. For more information
#'   about SomaScan identifiers and their formats, please see
#'   [SomaDataIO::SeqId].
#' @author Amanda Hiser
#' @examples
#' path2apt("GO:0019319")
#'
#' # Multiple pathways can be mapped at once
#' path2apt(c("M41209", "GO:2001256", "M1568"))
#'
#' # By default, pathway genes are matched to analytes by gene symbol. Matching
#' # by Entrez Gene ID instead can recover analytes whose gene was renamed
#' # (e.g. DDX58 -> RIGI)
#' sym <- path2apt("M5890")
#' ent <- path2apt("M5890", id_type = "EntrezGeneID")
#' setdiff(ent[[1L]], sym[[1L]])
#' @importFrom SomaDataIO getAnalyteInfo
#' @export
path2apt <- function(x,
                     col_meta_df = NULL,
                     id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                     verbose = interactive()) {

    id_type <- match.arg(id_type)

    stopifnot("`x` must be a character vector." = is.character(x))

    missing_ids <- setdiff(x, pathway_map$pathway_id)
    if ( length(missing_ids) > 0L ) {
        stop("Pathway ID(s) not found in `pathway_map`: ",
             paste(missing_ids, collapse = ", "), call. = FALSE)
    }

    if ( is.null(col_meta_df) ) {
        if ( verbose ) {
            message("`col_meta_df` not provided, using 11K annotations by default.")
        }
        col_meta_df <- SomaDataIO::getAnalyteInfo(example_data_11k)
    } else {
        stopifnot("`col_meta_df` must be a data frame or tibble." = inherits(col_meta_df, "data.frame"))
    }

    id_col  <- .match_col(col_meta_df, id_type)
    map_col <- ifelse(id_type == "EntrezGeneSymbol", "gene_symbol", "entrez_id")
    id_map  <- .make_id_map(col_meta_df, 
                            gene_col = id_col, 
                            apt_col = "AptName")

    path_df    <- pathway_map[pathway_map$pathway_id %in% x, c("pathway_id", map_col)]
    path_genes <- split(path_df[[map_col]],
                        factor(path_df$pathway_id, levels = unique(x)))

    apt_list <- lapply(path_genes, function(genes) {
        unique(id_map$AptName[id_map[[id_col]] %in% genes])
    })

    return(apt_list)
}
