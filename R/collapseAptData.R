#' Transform Aptamer-Centric Data to Gene-Centric Data
#'
#' Utility function to convert a `SeqId` or `AptName`-centric data frame, with 
#' one row per analyte, to a gene-centric data frame, with one row per gene. 
#' In cases where more than one analyte is associated with a given gene, the 
#' analytes will be collapsed into one row, with a delimiter to separate unique 
#' values.
#' 
#' @inheritParams apt2gene
#' @param df Data frame to be converted, most often the output of 
#'   [SomaDataIO::getAnalyteInfo()]. Must contain, at minimum, an "AptName"
#'   column, and the column specified by `id_type`.
#' @param apt_order Character. The name of a column in `df` containing 
#'   **numeric** ordering information for collapsed SomaScan identifiers. 
#'   By default, no ordering is used. Should only be specified when multiple 
#'   aptamers are associated with a given gene, and a particular order is 
#'   desired when the values are pasted together into a single, delimited 
#'   string.
#' @param sep Character used to separate collapsed SomaScan identifiers.
#'   Default is "|". If NULL, values will *not* be collapsed, and a character 
#'   vector will be returned in each column field when >1 SomaScan identifier 
#'   is associated with a given gene. See examples.
#' @returns A data frame with one row per unique gene.
#' @author Amanda Hiser, Alex Poole
#' @examples
#' df <- SomaDataIO::getAnalyteInfo(example_data_11k)
#' res <- collapseAptData(df = df)
#' 
#' # Pulling a couple example genes
#' res[res$EntrezGeneSymbol %in% c("NOTCH1", "APOE", "CAMP"), ]
#'
#' # Specifying a different delimiter
#' res_delim <- collapseAptData(df, sep = ", ")
#' res_delim[res_delim$EntrezGeneSymbol %in% c("NOTCH1", "APOE", "CAMP"), ]
#'
#' # Using a ranking metric to reorder the analytes
#' df$Analyte_Rank <- sample(1:nrow(df), nrow(df), replace = FALSE)
#' res_ord <- collapseAptData(df, apt_order = "Analyte_Rank")
#' res_ord[res_delim$EntrezGeneSymbol %in% c("NOTCH1", "APOE", "CAMP"), ]
#' 
#' #' Setting `sep = NULL` to create a list-column of SomaScan IDs
#' res_list <- collapseAptData(df, apt_order = "Analyte_Rank", sep = NULL)
#' res_list$AptName["APOE"]
#' @export
collapseAptData <- function(df, 
                            id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                            sep = "|", 
                            apt_order = NULL) {
    
    id_type <- match.arg(id_type)
    id_col  <- .match_col(df, id_type)

    # Create the id_map using the utility function
    id_map <- .make_id_map(col_meta_df = df, 
                           gene_col = id_col, 
                           apt_col = "AptName")
    
    # Add ordering info to id_map, if specified
    if ( !is.null(apt_order) ) {
        if ( !is.numeric(df[[apt_order]]) && !is.integer(df[[apt_order]]) ) {
            warning(paste("Ordering may produce unexpected results when the",
                          "`apt_order` column is not numeric."))
        }
        
        # Preserve order of original data frame
        ord_lookup <- setNames(df[[apt_order]], df[["AptName"]])
        id_map$id_order <- ord_lookup[id_map$AptName] # Re-order to match
        id_map <- id_map[order(id_map$id_order), ]
        id_map$id_order <- NULL # No longer needed in data frame
    }
    
    if ( !is.null(sep) ) {
        # Outputs a vector, will need to make a df with the results
        collapsed <- .collapse_by_gene(id_map, 
                                       apt_col = "AptName", 
                                       gene_col = id_col, 
                                       sep = sep)
        unique_gene_df <- data.frame(gene = names(collapsed),
                                     AptName = unname(collapsed))
        colnames(unique_gene_df)[1] <- id_col # Rename to match actual column name
    } else {
        # When sep == NULL, return list-columns
        split_df <- split(id_map, id_map[[id_col]])
        gene_list <- lapply(split_df, function(x) unique(x[["AptName"]]))
        
        unique_gene_df <- data.frame(gene = names(gene_list), 
                                     stringsAsFactors = FALSE)
        unique_gene_df$AptName <- gene_list
        colnames(unique_gene_df)[1] <- id_col
    }
    
    return(unique_gene_df)
}
