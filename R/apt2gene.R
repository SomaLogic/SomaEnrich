#' Convert SomaScan Analytes to Genes
#'
#' Transform a SomaScan analyte into the gene associated with the analyte's
#' target.
#'
#' @param x Character. A vector of SomaScan analyte identifiers in `AptName`
#'   (seq.1234.56) format. Blank ("") values will be coerced to
#'   NA. For more information about SomaScan identifiers and their formats,
#'   please see [SomaDataIO::SeqId].
#' @param col_meta_df Data frame of column metadata, likely created by
#'   [SomaDataIO::getAnalyteInfo()]. This object should contain SomaScan menu
#'   information that can be used to map between `AptNames` and Entrez gene
#'   identifiers. At minimum, it must contain columns named "AptName"
#'   (containing SomaScan identifiers in `AptName` format) and either
#'   "EntrezGeneSymbol" or "EntrezGeneID".
#'   
#'   If not specified, column metadata will be retrieved from `example_data_11k`
#'   by default.
#' @param id_type Character. Type of gene identifier to be used in the output.
#'   Options are "EntrezGeneSymbol" or "EntrezGeneID". Default is
#'   "EntrezGeneSymbol".
#' @param collapse Logical. Should values be collapsed to maintain the length of
#'   the original vector? Default is TRUE. Most useful when the results will be
#'   used as a data frame column. When FALSE, all genes associated with a given
#'   `AptName` will be returned as individual elements in the output vector. See
#'   examples.
#' @param verbose Logical. Should function messages be printed to the console?
#' @returns A character vector of genes associated with the values provided to
#'   `x`. `NA` is returned when a match is not found.
#' @author Amanda Hiser, Alex Poole
#' @seealso [SomaDataIO::getAnalyteInfo()]
#' @examples
#' anno <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
#' apts <- withr::with_seed(123,
#'                          sample(SomaDataIO::getAnalytes(SomaDataIO::example_data), 10))
#' apts <- c(apts, "seq.10367.62") # Adding IL12
#' apt2gene(x = apts, col_meta_df = anno)
#'
#' # Adding results to a data frame
#' df <- withr::with_seed(321, data.frame(AptName = apts,
#'                                        Value = rnorm(11)))
#' df$Gene <- apt2gene(x = df$AptName, col_meta_df = anno)
#' df
#'
#' # Using `collapse = FALSE` to separate all returned values
#' apt2gene(x = df$AptName, col_meta_df = anno, collapse = FALSE)
#' @importFrom SomaDataIO is.AptName getAnalyteInfo
#' @export
apt2gene <- function(x,
                     col_meta_df = NULL,
                     id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                     collapse = TRUE,
                     verbose = interactive()) {

    id_type <- match.arg(id_type)

    if ( inherits(x, "list") ) {
        stop("`x` must be a vector, not a list", call. = FALSE)
    }

    if ( !all(is.AptName(x)) ) {
        stop("At least some values in `x` are not 'AptNames': ",
             paste(x[!is.AptName(x)], collapse = ", "),
             call. = FALSE
        )
    }

    if ( is.null(col_meta_df) ) {
        if ( verbose ) {
            message("`col_meta_df` not provided, using 11K annotations by default.")
        }
        col_meta_df <- SomaDataIO::getAnalyteInfo(example_data_11k)
    }

    stopifnot(
        "`col_meta_df` must be a data frame, tibble, or similar." = inherits(col_meta_df, "data.frame")
    )

    id_col <- .match_col(col_meta_df, id_type)
    genes <- col_meta_df[match(x, col_meta_df$AptName), ]
    genes <- genes[[id_col]]

    if ( collapse ) {
        genes
    } else {
        strsplit(genes, "\\||\\,\\s|\\s+") |>
            unlist() |>
            trimws(which = "both")
    }
}
