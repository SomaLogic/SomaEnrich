#' Select A Single Measurement Per Feature
#'
#' @inheritParams somaPrGSEA
#' @inheritParams .select_best
#' @param verbose Logical. Should informative progress messages be printed to 
#'   the console?
#' @returns A list containing two elements:
#'  \item{results}{The filtered input vector, with non-prioritized features removed.}
#'  \item{removed}{A named vector containing elements removed from the input 
#'                 `ranks` vector.}
#' @importFrom tidyr separate_rows
#' @importFrom SomaDataIO getAnalyteInfo
#' @importFrom stats setNames
#' @keywords internal
.resolve_many_to_one <- function(ranks, 
                                 resolve_method = c("abs", "min", "max", "rank"),
                                 verbose = interactive()) {
    
    resolve_method <- match.arg(resolve_method)
    
    top_ranked <- .select_best(ranks, resolve_method = resolve_method)
    keep_idx <- unname(unlist(top_ranked$keep_list))
    keep_idx <- sort(keep_idx)
    
    # Retrieves values of non-highest ranked list element(s)
    rmv <- top_ranked$ranks[-keep_idx]
    
    # Filters id_list to only features w/ multimapping present
    multi <- top_ranked$id_list[lengths(top_ranked$id_list) > 1]
    
    if ( verbose ) {
        n_multi <- length(multi)
        n_rmv <- length(rmv)
        .todo("There are {.val {n_multi}} feature(s) with multiple elements in the ranked input vector.")
        .done("A total of {.val {n_rmv}} lower-ranked element(s) were dropped.")
    }
    
    final_ranks <- top_ranked$ranks[keep_idx]

    list(results = final_ranks,
         removed = rmv)
}

