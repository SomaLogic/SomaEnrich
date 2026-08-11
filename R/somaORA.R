#' Perform Overrepresentation Analysis for SomaScan
#'
#' @description
#'   Implements the algorithm from [fgsea::fora()] to perform
#'   overrepresentation analysis (ORA) via hypergeometric test, tailored
#'   specifically to proteomic SomaScan data.
#'
#' @section Additional Details:
#'   Overrepresentation analysis is a statistical method that determines
#'   whether sets of features, often genes, are present in a data set
#'   more than would be expected by chance. As an example, imagine a SomaScan
#'   experiment in which 200 analytes are determined to be differentially
#'   expressed in a group of breast cancer samples. Using ORA, a
#'   researcher can investigate whether any pathways, like cellular
#'   proliferation, are more highly represented in those 200 analytes than
#'   would be expected by chance.
#'
#' @section Feature Selection: 
#'   The primary input to ORA is a set of biological features of interest.
#'   Selecting these features typically involves defining a specific subset,
#'   often differentially expressed genes, based on previously defined
#'   thresholds or cutoffs. For example, one could use the combined statistical
#'   significance and log fold change (logFC) cutoffs of p < 0.05 and
#'   logFC > 2, respectively, to identify significantly differentially
#'   expressed genes. This list would then be used as input for ORA. Please
#'   note that the most appropriate cutoff(s) to use will likely vary by
#'   experiment; this is only an example.
#'
#' @section Enrichment Universe:
#'   The `universe` (sometimes referred to as "background") serves as a point
#'   of reference against which the genes or analytes of interest (provided to
#'   `features`) can be compared. Setting the correct background universe is
#'   _crucial_ for accurate ORA results. The universe must encompass all
#'   features in the analysis results (typically differential expression
#'   analysis). If analyte filtering was performed prior to differential
#'   expression analysis, the filtered analytes should not be used to set the
#'   ORA universe.
#'
#'   If `features` contains elements that are not in `universe`, those
#'   elements will be dropped from `features`. Additionally, `universe`
#'   elements that are not found in any gene set will be dropped. 
#'
#' @section Handling Heterodimers:
#'   When a SOMAmer Reagent targets a protein heterodimer, the genes associated
#'   with that heterodimer are collapsed using a delimiter character
#'   (ex. "CGA|CGB3|CGB7") for concise display in the SomaScan menu. However,
#'   when performing ORA, these collapsed and delimited strings must be split
#'   into individual gene symbols, so each symbol can be matched against
#'   gene-based databases like GO or MSigDB. In SomaEnrich, the `resource` or
#'   `cust_paths` arguments direct annotation of the features of interest to
#'   known pathways or gene sets. If the elements of `features` do not match
#'   elements in `resource`/`cust_paths`, the non-matching elements in
#'   `features`/`universe` will be dropped. Consequently, to preserve these
#'   elements, the default behavior is to split delimited strings.
#'
#' @inheritParams somaPrGSEA
#' @inheritSection somaPrGSEA Min/Max Set Sizes
#' @param features Character. A vector containing features of interest, often
#'   differentially expressed genes or SomaScan analytes. Vector elements must
#'   be either Entrez gene symbols/identifiers (unless custom gene sets are 
#'   used, see `cust_paths`) or SomaScan identifiers in R-compatible 
#'   `AptName` format (e.g. `seq.1234.56`). If using `AptNames`, you must set 
#'   `use_aptnames = TRUE`.
#' @param id_type Character. Type of gene identifier used in `features`.
#'   Options are "EntrezGeneSymbol" or "EntrezGeneID". Default is 
#'   "EntrezGeneSymbol". Ignored if `use_aptnames = TRUE`.
#' @param split_ids Logical. If the identifiers used in the `features` and
#'   `universe` vectors contain heterodimers (elements separated by a
#'   delimiter, like `"GENE1|GENE2"`), should the heterodimers be split into
#'   individual components, like `c("GENE1", "GENE2")`? If TRUE (the default),
#'   this may introduce duplicate identifiers into `features`/`universe`.
#'   Applies to _both_ `features` and `universe`. See `Details` and examples.
#' @param unique_features Logical. If present, should duplicates in
#'   `features` be removed? Default is FALSE. If TRUE, the resultant
#'   vector will only contain unique elements. The duplicated features are
#'   often the result of setting `split_ids = TRUE` when heterodimers
#'   are present in the `features` vector. See examples.
#' @param universe Character. A vector of unique gene symbols or `AptNames`.
#'   This vector represents the overall "universe" from which the elements in
#'   `features` were identified. The `universe` should encompass _all_
#'   elements in the differential expression analysis from which the set of
#'   over-represented targets (`features`) was chosen. See `Details` for more
#'   information.
#'   If using `AptNames`, you must set `use_aptnames = TRUE`.
#' @returns A `data.frame` object containing the results of ORA, where the
#'   results in each row correspond to a single tested pathway. The data frame
#'   contains the following columns:
#' \item{resource_code}{Abbreviated character string representing the name 
#'                      of the original gene set resource.}
#' \item{pathway_id}{Pathway identifier/accession number from the original
#'                   gene set or pathway resource.}
#' \item{pathway}{Name of the pathway.}
#' \item{pval}{An enrichment p-value from hypergeometric test.}
#' \item{padj}{The BH-adjusted p-value.}
#' \item{foldEnrichment}{The degree of enrichment, relative to the background.}
#' \item{overlap}{The size of the overlap between the pathway and features of interest.}
#' \item{size}{The size of the pathway.}
#' \item{overlapFeatures}{A vector containing the overlapping features (see 'overlap' column).}
#'
#' By default, the data frame is sorted by the BH-adjusted p-value and fold
#' enrichment.
#' @author Amanda Hiser
#' @examples
#' # Genes mapping to differentially expressed analytes
#' deg <- t_tests[t_tests$fc > 1 & t_tests$fdr < 0.05, ]$EntrezGeneSymbol
#' 
#' # Input universe = all genes in statistical test results
#' uni <- t_tests$EntrezGeneSymbol
#'
#' # To properly check for enrichment, values in the input vector must match 
#' # values in the gene sets/pathway resource. 
#' # SomaScan contains heterodimers that won't match non SomaScan-derived gene sets
#' head(deg)
#' deg[grepl("\\|", deg)]
#'
#' # Splitting these values can create duplicates in the input vector
#' split_deg <- strsplit(deg, "\\|") |> unlist()
#' split_deg[duplicated(split_deg)]
 #' 
#' # The following will run ORA with GO BP gene sets (the default),
#' # split heterodimers into individual components (the default),
#' # and remove duplicates created from splitting
#' res <- somaORA(features = deg,
#'                universe = uni,
#'                unique_features = TRUE)
#'
#' # ORA can be performed using SomaScan AptNames, instead of genes
#' dea <- head(t_tests$AptName, 50)
#' apt_uni <- t_tests$AptName
#'
#' # A progress bar will appear to track AptName conversion
#' meta <- SomaDataIO::getAnalyteInfo(example_data_11k)
#' res_apt <- somaORA(features = dea,
#'                    use_aptnames = TRUE,
#'                    universe = apt_uni,
#'                    col_meta_df = meta)
#' @importFrom fgsea fora
#' @importFrom SomaDataIO is.AptName
#' @export
somaORA <- function(features,
                    universe,
                    resource = c("bp", "mf", "h", "c1", "c3",
                                 "c4", "c6", "c7", "c8"),
                    split_ids = TRUE,
                    unique_features = TRUE,
                    id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                    use_aptnames = FALSE,
                    col_meta_df = NULL,
                    cust_paths = NULL,
                    min_feats = 5L,
                    max_feats = 500L,
                    verbose = interactive()) {

    resource <- match.arg(resource)
    resource <- ifelse(!is.null(cust_paths), "custom", resource)
    id_type  <- match.arg(id_type)
    id_col <- switch(id_type, 
                     EntrezGeneID = "entrez_id", 
                     EntrezGeneSymbol = "gene_symbol")
    
    # Making sure universe is cleaned up
    universe <- unique(universe)
    universe <- universe[universe != ""]
    universe <- universe[!is.na(universe)]

    if ( length(features) == 0L ) {
        stop("`features` must contain at least one element.", call. = FALSE)
    }
    
    if ( any(SomaDataIO::is.AptName(features)) && any(SomaDataIO::is.AptName(universe)) && !use_aptnames & is.null(cust_paths) ) {
        stop(
            "`features` and `universe` appear to be SomaScan identifiers,",
            "did you set `use_aptnames = TRUE`?",
            call. = FALSE
        )
    }

    if ( split_ids ) {
        features <- .splitIDs(features)
        universe <- .splitIDs(universe)
        
        if ( unique_features ) {
            features <- unique(features)
        }
    }
    
    if ( any(duplicated(universe)) ) {
        if ( verbose ) {
            writeLines(.rule("Pre-Processing", line_col = "magenta"))
            .oops("Duplicate elements found in {.code universe}. ")
            .todo("These elements will be removed to create a unique vector.")
        }
        # Can't have duplicate values in universe
        universe <- unique(universe)
        universe <- universe[universe != ""]
    }

    stopifnot("All elements in `features` must be found in `universe`." = all(features %in% universe))

    if ( !is.null(cust_paths) ) {
        stopifnot("`cust_paths` must be a list of named vectors." = inherits(cust_paths, "list"))
        if ( use_aptnames ) {
            warning("`use_aptnames` is ignored when `cust_paths` is provided. ",
                    "Custom pathways must already use the same identifiers as `features`.",
                    call. = FALSE)
        }
        sel_paths <- cust_paths
    } else {
        # Filter pathway map to retain resource specified by user
        sel_paths <- pathway_map[pathway_map$group_code == resource, ]
        
        # Add pathway name from addl pathway info lookup table
        name_map <- setNames(path_name_lookup$pathway_name, path_name_lookup$pathway_id)
        sel_paths$pathway_name <- unname(name_map[sel_paths$pathway_id])

        # Convert to list of named vectors for input into ORA
        sel_paths <- df2list(sel_paths, "pathway_name", id_col, clean_names = FALSE)
        
        if ( use_aptnames ) {
            stopifnot(
                "All features must be AptNames when `use_aptnames = TRUE`." = all(SomaDataIO::is.AptName(features)),
                "`col_meta_df` must be provided when `use_aptnames = TRUE`." = !is.null(col_meta_df)
            )

            # Convert gene paths to AptNames
            sel_paths <- genePath2aptPath(sel_paths,
                                          col_meta_df = col_meta_df,
                                          id_type = "EntrezGeneSymbol",
                                          verbose = verbose)
        }
    }

    if ( !any(features %in% unique(unlist(sel_paths))) ) {
        stop("No elements in `features` were found in the selected gene sets.",
             call. = FALSE)
    }

    if ( verbose ) {
        writeLines(.rule("Performing Enrichment Analysis", line_col = "cyan"))
        .todo("Running ORA...")
    }

    results <- fgsea::fora(pathways = sel_paths,
                           genes    = features,
                           universe = universe,
                           minSize  = min_feats,
                           maxSize  = max_feats)

    results <- results[order(results$padj, -abs(results$foldEnrichment), decreasing = FALSE), ]

    # Change default column name to be more general
    apt_idx <- which(colnames(results) == "overlapGenes")
    colnames(results)[apt_idx] <- "overlapFeatures"

    # Add pathway ID
    id_map  <- setNames(path_name_lookup$pathway_id, path_name_lookup$pathway_name)
    new_col <- data.frame(resource_code = resource,
                          pathway_id    = unname(id_map[results$pathway]))
    
    if ( verbose ) {
        .done("Done!")
    }
    
    return(cbind(new_col, results))
}
