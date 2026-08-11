#' Convert SomaScan ID-Based Feature Set into Gene Set
#'
#' @description
#' Utility function to convert a set of SomaScan identifiers (i.e. `AptNames`)
#' into genes, while also filtering to avoid potential gene
#' duplication produced by multiple SomaScan analytes mapping to the same gene.
#' 
#' Mapping between genes and SomaScan analytes is not always 1:1, and it's
#' possible for multiple analytes to map to the same gene. For example,
#' `13699-6`, `4330-4`, and `8468-19` all map to prostate-specific antigen
#' (PSA). If all of these analytes were included in a feature set, PSA would be
#' represented 3 times in that set. However, when performing over-representation
#' analysis, pathway analysis, or GSEA, it's crucial to avoid instances where a
#' single feature is represented multiple times in the same feature set. This
#' can lead to artificial over-representation when the duplicated feature is
#' counted multiple times.
#' 
#' After converting SomaScan identifiers to their associated gene, this
#' function will ensure that each gene is represented only _once_ in the gene
#' set or pathway. However, because 1:many mapping may be present between
#' identifiers, the pathway returned by `aptPath2genePath()` may be longer than
#' the input aptamer list.
#' Because of this feature, it is recommended to use this function instead of
#' `apt2gene()` when converting entire gene sets or pathways. See the examples
#' for more details.
#'    
#' @param path_list Named list containing character vectors of SomaScan
#'   `AptNames`.
#' @param col_meta_df Data frame of column metadata, likely created by
#'   [SomaDataIO::getAnalyteInfo()]. This object should contain SomaScan menu
#'   information that can be used to map between `AptNames` and Entrez gene
#'   identifiers. At minimum, it must contain columns named "AptName"
#'   (containing SomaScan identifiers in `AptName` format) and either
#'   "EntrezGeneSymbol" or "EntrezGeneID".
#' @param id_type Character. Desired gene identifier to be used in the output.
#'   Options include "EntrezGeneSymbol" or "EntrezGeneID". Default is
#'   "EntrezGeneSymbol".
#' @returns The input `pathway` vector, modified to use gene identifiers
#'   instead of `AptNames`. `AptNames` that do not have an associated gene in
#'   the SomaScan menu (determined using the column metadata in the `soma_adat`
#'   object provided to `adat =`) will be dropped.
#' @param verbose Logical. Should progress messages be printed to the console?
#' @author Amanda Hiser
#' @examples
#' # Gene symbols are returned by default
#' aptPath2genePath(path_list   = list(Example = ex_apt_pathway),
#'                  col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k))
#' 
#' # When a SomaScan assay analyte is associated with more than one gene,
#' # all genes will be returned
#' aptPath2genePath(path_list = list(Single_Analyte = "seq.16927.9"),
#'                  col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k))
#' @importFrom SomaDataIO is.AptName
#' @export
aptPath2genePath <- function(path_list, 
                             col_meta_df,
                             id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                             verbose = interactive()) {

  id_type <- match.arg(id_type)
  id_col  <- .match_col(col_meta_df, id_type)

  # Won't catch logical(0)
  stopifnot("All SomaScan identifiers should be in 'AptName' format." = all(SomaDataIO::is.AptName(unlist(path_list))))

  # Unique vector of all features in provided pathway, converted to AptNames
  all_apts <- unique(unname(unlist(path_list)))

  # Retrieving annotations to map analytes -> genes
  assay_anno <- .make_id_map(col_meta_df, gene_col = id_col)

  if ( !all(all_apts %in% assay_anno$AptName) && verbose ) {
    message(paste0(sum(!(all_apts %in% assay_anno$AptName)),
                   " analytes in 'path_list' do not have an associated gene",
                   " identifier in 'col_meta_df'.\n",
                   "These analytes will be dropped."))
  }

  result <- lapply(path_list, function(x) {
    path_anno <- assay_anno[assay_anno$AptName %in% x, ]

    # Map AptNames to user-specified gene ID type
    conv_path <- unique(path_anno[[id_col]])

    # Final cleanup of converted values
    conv_path <- conv_path[!is.na(conv_path)]
  })

  return(result)
}
