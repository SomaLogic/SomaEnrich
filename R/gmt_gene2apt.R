#' Convert Gene-Based GMT File to AptNames
#' 
#' Reads in a standard GMT file (using gene identifiers) and converts genes
#' to SomaScan `AptNames`. 
#' 
#' @inheritParams genePath2aptPath
#' @param gmt_gene_file Path to the `.gmt` file to be converted.
#' @param write Logical. Should the converted data set be written out to 
#'    a new GMT file? Default is FALSE.
#' @param gmt_apt_file The desired file name for the converted file. Required 
#'    if `write = TRUE`. The file will be written to the current directory, 
#'    unless a path to the final file is provided.
#' @returns A named list of gene sets/pathways, converted to SomaScan 
#'    `AptNames`.
#' @examples
#' \dontrun{
#' example_gmt <- system.file("extdata", 
#'                            "msigdb_h_androgen_response.gmt", 
#'                             package = "SomaEnrich")
#' gmt_gene2apt(example_gmt, 
#'              col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k), 
#'              gmt_apt_file = "msigdb_h_androgen_response_apt.gmt")
#' }
#' @importFrom fgsea gmtPathways writeGmtPathways
#' @export
gmt_gene2apt <- function(gmt_gene_file, 
                         col_meta_df,
                         id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                         write = FALSE,
                         gmt_apt_file = NULL,
                         verbose = interactive()) {
    
    id_type <- match.arg(id_type)
    
    x <- fgsea::gmtPathways(gmt_gene_file)
    
    converted <- genePath2aptPath(x, id_type  = id_type, 
                                  col_meta_df = col_meta_df,
                                  verbose     = verbose)
    
    if ( write ) {
      
      if ( is.null(gmt_apt_file) ) {
        stop("Please provide a filename to `gmt_apt_file` when `write = TRUE`.")
      }
      
      if ( verbose ) {
        message(paste("Writing converted GMT file to", gmt_apt_file)) 
      }
      fgsea::writeGmtPathways(converted, gmt_apt_file)
    }
    
    return(converted)
}
