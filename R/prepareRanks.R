#' Prepare Ranking Metrics for SomaScan-Based Enrichment Analyses
#' 
#' @description
#'   Resolves instances of identifier duplication in SomaScan datasets,
#'   typically resulting from many-to-one relationships between SOMAmer
#'   Reagents and their protein targets. This duplication is often encountered
#'   when preparing a statistical ranking metric for pre-ranked GSEA, or
#'   (more generally) converting SomaScan analyte measurements/statistics to
#'   the gene level.
#'   
#'   This function will ensure that the provided data only has one numeric
#'   value per unique feature (typically genes), and in cases of duplication,
#'   a single value will be selected to represent each gene.
#' 
#' @section SOMAmer Reagent and Target Relationships:
#'   Non-1:1 mapping between SOMAmer Reagents and their targets is a known
#'   characteristic of the SomaScan menu. This can lead to a single gene being
#'   associated with multiple statistical metric results. For example,
#'   prostate-specific antigen (PSA) and its associated gene, KLK3, are
#'   targeted by 4 SOMAmer Reagents in the 7K SomaScan menu. 
#'   This means that, _at the gene level_, KLK3 has 4 RFU values in any 7K
#'   SomaScan data set.
#' 
#'   In pre-ranked GSEA, however, values in the input ranking data must be
#'   unique. Therefore, the input data cannot contain 4 measurements for PSA -
#'   only one must be chosen to be representative of the KLK3 gene. There are
#'   many acceptable selection methods, and the most appropriate one will
#'   depend on your question and/or experimental design. See the 
#'   analyte-to-gene mapping vignette (`browseVignettes(package = "SomaEnrich")`)
#'   for more details.
#'   
#' @section Ties in the Ranking Metric:
#'   When a SOMAmer Reagent targets a protein heterodimer, the genes associated
#'   with that heterodimer are collapsed using a delimiter character
#'   (ex. "CGA|CGB3|CGB7") for concise display in the SomaScan menu. However,
#'   when performing GSEA, the collapsed and delimited gene identifier must
#'   be split into individual gene symbols, so each symbol can
#'   be matched against gene-based databases like GO or MSigDB. Splitting a
#'   single analyte into multiple genes requires duplicating the ranking metric
#'   value across all resultant genes (ex. `c("GENE1|GENE2" = 1)` becomes
#'   `c("GENE1" = 1, "GENE2" = 1`). This can introduce ties into the ranked
#'   input vector; however, this is _expected_ behavior. The tied values
#'   indicate that there is insufficient information to rank the genes that
#'   were co-measured as part of the same protein complex. See the 
#'   analyte-to-gene mapping vignette (`browseVignettes(package = "SomaEnrich")`)
#'   for more details.
#'   
#' @inheritParams .resolve_many_to_one
#' @inheritParams .select_best
#' @param stats Vector of a numeric statistic. This vector will
#'   be used to rank features, according to the value of the statistic. Any
#'   numeric metric may be used (t-statistic, p-value, fold change, etc.).
#'   When `resolve_method = "rank"`, the most important or highly valued
#'   features should be sorted _first_. This may require values to be
#'   ordered in ascending or descending order, depending on the metric used.
#'   For example, p-values should be sorted in ascending order, with smaller
#'   values indicating higher importance, while fold-change values could be
#'   sorted in descending order, with larger positive values indicating higher
#'   importance.
#' @param features Character. Identifiers that correspond to the statistical
#'   metric in `stats`, often SomaScan `AptNames` or genes.
#' @param split_ids Logical. If the identifiers used for the statistical
#'   metric contain heterodimers (features separated by a delimiter, like
#'   `"GENE1|GENE2"`), should the heterodimers be split into individual
#'   components, like `c("GENE1", "GENE2")`? If TRUE (the default), the final
#'   ranking metric will contain duplicated numeric values, one for each
#'   unique member of the heterodimer. See `Details` and examples.
#' @param resolve_multimapping Logical. Should many-to-one relationships between
#'   SOMAmer Reagents and genes (i.e. multiple SOMAmer Reagents mapping to the
#'   same gene) in the statistical ranking metric be resolved? Default is TRUE.
#'   Only one entry per gene will be retained, based on the `resolve_method`
#'   selection criteria. See `Details`.
#' @param verbose Logical. Should progress messages be printed to the console?
#'   By default, messages will only be shown in interactive R sessions. Set
#'   to TRUE to show all messages, and FALSE to silence messages.
#' @returns A named numeric vector of statistics, with unique element
#'   names. This vector is suitable for input into `somaPrGSEA()`, or other
#'   gene-based analytical methods.
#' @examples
#' # In some cases, multiple SomaScan analytes map to the same gene ID
#' example_apt_vec <- c("seq.13699.6", "seq.21232.39", "seq.4330.4", "seq.8468.19")
#' apt2gene(example_apt_vec, SomaDataIO::getAnalyteInfo(example_data_11k))
#' 
#' # All KLK3 analytes are present in the t-tests results data, each with
#' # distinct t-test statistics and significance values
#' t_tests[t_tests$EntrezGeneSymbol == "KLK3", ]
#' 
#' # Pre-ranked GSEA input must contain *one* unique ranking value per gene -
#' # using the following as input will be problematic
#' ranks <- t_tests$t_stat
#' names(ranks) <- t_tests$EntrezGeneSymbol
#' head(ranks, n = 10) # LEP appears multiple times
#' 
#' # Identifiers used for the input vector must match those in the gene sets.
#' # SomaScan contains heterodimers that won't match non SomaScan-derived
#' # gene sets
#' ranks[grepl("\\|", names(ranks))][1:5]
#' 
#' # prepareRanks() solves both of these problems.
#' # The default resolve_method = "abs" selects the analyte with the largest
#' # absolute value per gene. No pre-sorting required.
#' unique_ranks <- prepareRanks(stats = t_tests$t_stat,
#'                              features = t_tests$EntrezGeneSymbol)
#' 
#' head(unique_ranks)
#' any(duplicated(names(unique_ranks)))
#' 
#' # Caveat: splitting up heterodimers can introduce ties.
#' # All genes split from the original heterodimer will have the same value
#' ranks["CGA|FSHB"]
#' unique_ranks[c("CGA", "FSHB")]
#' @export
prepareRanks <- function(stats,
                         features,
                         split_ids = TRUE,
                         resolve_multimapping = TRUE,
                         resolve_method = c("abs", "min", "max", "rank"),
                         verbose = interactive()) {

    resolve_method <- match.arg(resolve_method)

    if ( length(stats) != length(features) ) {
        stop("`stats` and `features` must have the same length.", call. = FALSE)
    }

    input_vec <- stats
    names(input_vec) <- features

    # Capture sort direction of original stats before any filtering.
    # Used later to re-sort output when resolve_method = "rank"
    desc_input <- is.unsorted(stats, na.rm = TRUE)

    if ( anyNA(input_vec) ) {
        n_na <- sum(is.na(input_vec))
        .warn("{.val {n_na}} NA value{?s} detected in `stats`. They may be dropped if `resolve_multimapping = TRUE`.")
    }
    
    if ( split_ids ) {
        # First check to see if heterodimers/complexes are present
        if ( any(grepl("\\|", features)) ) {
            if ( verbose ) {
                writeLines(.rule("Pre-Processing Inputs", line_col = "magenta"))
                .oops("Heterodimers are present in the ranking vector names.")
                .todo("They will be split into individual elements.")
            }
            n_original <- length(input_vec)
            input_vec <- .splitIDs(input_vec)
            n_split <- length(input_vec)
            if ( verbose ) {
                .done("The ranking vector was expanded from {.val {n_original}} to {.val {n_split}} elements.\n")
            }
        }
    }
    
    # Return early if no cases of non-1:1 mapping are present.
    # The heterodimers have been split and multimapping resolution doesn't have
    # anything to do
    if ( !any(duplicated(names(input_vec))) ) {
        return(input_vec)
    }
    
    if ( !resolve_multimapping ) {
        stop("Duplicate values found in the ranking vector names. ",
             "Please remove duplicates before proceeding, or set ",
             "'resolve_multimapping = TRUE'")
    }
    
    if ( verbose ) {
        writeLines(.rule("Resolving Non-Unique Mapping", line_col = "blue"))
        .oops("Duplicate value(s) found in the ranking vector names.")
        .todo("Filtering will be performed to retain a single value.")
    }

    result <- .resolve_many_to_one(input_vec,
                                   resolve_method = resolve_method,
                                   verbose = verbose)$results

    result
}
