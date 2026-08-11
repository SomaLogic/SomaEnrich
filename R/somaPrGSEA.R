#' Perform Pre-Ranked GSEA for SomaScan
#'
#' @description
#'   Wrapper function for [fgsea::fgsea()] to perform pre-ranked
#'   gene set enrichment analysis (GSEA), tailored specifically to proteomic
#'   SomaScan data.
#'
#' @section Additional Details:
#'   Gene Set Enrichment Analysis (GSEA) is a rank-based statistical method that
#'   tests whether predefined groups of biologically related entities are
#'   over-represented at the top or bottom of a ranked list of features.
#'   Unlike many other methods of biological interpretation (like ORA), GSEA
#'   does not rely on arbitrary significance cutoffs.
#'   In proteomics experiments, GSEA is used to detect coordinated changes
#'   across groups of functionally related proteins/genes, and provides a source
#'   for functional biological interpretation of statistical test results.
#'
#' @section Min/Max Set Sizes:
#'   The statistical power of GSEA/ORA to detect enrichment depends on the size
#'   of the gene set being tested, and GSEA normalization is not accurate
#'   for very small feature sets. The minimum allowable gene set is defined by
#'   the `min_feats` argument; any sets smaller or larger than `min_feats` or
#'   `max_feats` will be discarded _before_ performing GSEA/ORA. Please note
#'   that this threshold is applied _after_ filtering to remove features from
#'   each gene set that are not also present in the `ranks`, i.e. taking the
#'   intersect of the `ranks` and gene lists.

#'   The default `min_feats` threshold used here differs from other enrichment
#'   analysis tools. The 15-gene minimum suggested by the Broad Institute's
#'   [GSEA software](https://www.gsea-msigdb.org/gsea/index.jsp) was established
#'   in a whole-genome context, where the previously described intersection step
#'   removes few genes. However, in the reduced feature space of SomaScan,
#'   more genes are lost from taking this intersection. In some cases, this can
#'   greatly reduce the size of available gene sets. To account for this loss, a
#'   default minimum threshold of 5 is used. Consider final gene set size
#'   accordingly when reviewing your enrichment analysis results.
#'
#' @inheritParams prepareRanks
#' @inheritParams .select_best
#' @inheritSection prepareRanks SOMAmer Reagent and Target Relationships
#' @param ranks Either the output of `prepareRanks()`, or a named vector of a
#'   numeric ranking statistic. This vector will be used to rank features,
#'   according to the value of the statistic. Any numeric metric may be used
#'   (t-statistic, -log10(p), fold change, etc.). The vector _must_ be
#'   named using Entrez gene symbols/identifiers _or_ SomaScan identifiers in
#'   R-compatible `AptName` format (e.g. `seq.1234.56`). If using `AptNames`,
#'   set `use_aptnames = TRUE`.
#' @param id_type Character. Type of gene identifier used in `ranks`.
#'   Options are "EntrezGeneSymbol" or "EntrezGeneID". Default is
#'   "EntrezGeneSymbol". Ignored if `use_aptnames = TRUE`.
#' @param resource Character. String specifying the gene set database resource
#'   to use. Choose one of:
#'   \describe{
#'     \item{"bp"}{GO Biological Process}
#'     \item{"mf"}{GO Molecular Function}
#'     \item{"h"}{MSigDB Hallmark}
#'     \item{"c1"}{MSigDB Positional}
#'     \item{"c2"}{MSigDB Curated (note: does not contain KEGG or BioCarta)}
#'     \item{"c3"}{MSigDB Regulatory Target}
#'     \item{"c4"}{MSigDB Computational (note: does not contain CM submodule)}
#'     \item{"c6"}{MSigDB Oncogenic Signature}
#'     \item{"c7"}{MSigDB Immunologic Signature}
#'     \item{"c8"}{MSigDB Cell Type Signature}
#'   }
#'   The default is `bp`. For more advanced filtering, or to use gene sets not
#'   provided here, see the `cust_paths` argument. Note that `resource` is
#'   ignored when a data set is supplied to `cust_paths`.
#' @param cust_paths Optional. A named list containing character vectors of
#'   gene sets or functional groups of interest. If provided, `cust_paths` will
#'   override `resource`. The name of each list element should be the name of
#'   the pathway/gene set. The identifier used in each vector should be
#'   the same type as the names of `ranks`. This object is still subject to the
#'   minimum/maximum feature cutoffs defined by `min_feats` and `max_feats`.
#'
#'   This argument is meant to be used when the desired biological network is
#'   not present in the `resource` options.
#' @param use_aptnames Logical. Should analysis be performed in aptamer-space,
#'   i.e. using SomaScan `AptNames` instead of genes? Default is FALSE. If TRUE,
#'   the input ranking data must contain `AptNames` instead of gene
#'   identifiers. If using `cust_paths` when `use_aptnames = TRUE`, custom
#'   gene set members must also be `AptNames`.
#' @param col_meta_df Optional. Data frame of column metadata, likely created
#'   by calling [SomaDataIO::getAnalyteInfo()] on the `soma_adat` used for
#'   differential expression analysis. Required when
#'   `use_aptnames = TRUE`. This object will be used to map between SomaScan
#'   `AptNames` and gene symbols.
#' @param min_feats Numeric. Minimum number of features required for a gene set
#'   to be retained. Default is 5. Sets smaller than this value will be
#'   discarded. If `NULL`, no minimum is used. See `Details` for more
#'   information.
#' @param max_feats Numeric. Maximum number of features required for a gene set
#'   to be retained. Default is 500. See `Details` for more information.
#' @param n_perm Numeric. Number of permutations to perform. Default is 1000.
#' @param seed Numeric. Seed to set for GSEA calculations. Because the GSEA
#'   algorithm used here is non-deterministic, a seed must be set for results
#'   to be reproducible. Default is 42.
#' @param ... Optional arguments passed to [fgseaMultilevel()].
#' @return A list with two elements:
#' \describe{
#'   \item{results}{A `data.frame` object containing the results of GSEA,
#'                  where the results in each row correspond to a single tested
#'                  gene set.}
#'   \item{final_ranks}{The final ranks used as input for GSEA. If
#'                      `resolve_multimapping = TRUE` and `split_ids = TRUE`, the
#'                      ranks used as input were modified to split protein
#'                      heterodimers into unique gene members and resolve
#'                      non-1:1 mapping. Therefore, this vector may differ
#'                      from what was provided as input to `somaPrGSEA()`.}
#' }
#' 
#' The `results` data frame contains the following columns:
#' \item{resource_code}{Abbreviated character string representing the name
#'                      of the original gene set resource.}
#' \item{pathway_id}{Pathway identifier/accession number from the original
#'                   gene set or pathway resource.}
#' \item{pathway}{Name of the pathway.}
#' \item{pval}{Enrichment p-value.}
#' \item{padj}{FDR-adjusted p-value.}
#' \item{log2err}{The expected error for the standard deviation of the p-value
#'                logarithm.}
#' \item{ES}{The enrichment score calculated via the 
#'           [Broad GSEA implementation](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Enrichment_Score_(ES)).
#'           The ES represents the degree to which a set `S` is over-represented
#'           at the top or bottom of a ranked list of genes (`L`). The ES is
#'           calculated by walking down the list `L`, increasing a running-sum
#'           statistic whenever a gene in `L` is found in set `S`, and decreasing
#'           it each time a gene is not found in `S`. The final ES is the maximum
#'           deviation from zero encountered in the random walk. A positive ES
#'           indicates that `S` is enriched toward the top of the ranked
#'           list, while a negative ES indicates that `S` is enriched at the
#'           bottom of the ranked list.}
#' \item{NES}{Enrichment score normalized to account for the size of each gene
#'            set, calculated by taking the ES and dividing it by the mean ES
#'            for all data set permutations. The NES is the primary statistic 
#'            to use for examining GSEA results.}
#' \item{leadingEdge}{Vector with indexes of leading edge features that drive the
#'                    enrichment, see [Running a Leading Edge Analysis](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading).
#'                    Leading edge features in the set `S` appear in the ranked
#'                    list `L` at, or before, the point where the running sum
#'                    reaches its maximum deviation from zero. Essentially,
#'                    the leading edge subset is the core of a feature set that
#'                    accounts for the enrichment signal.}
#' \item{starting_set_size}{Original size of the pathway/gene set, without modifications.}
#' \item{final_set_size}{Number of SomaScan analytes mapping to
#'                       features in the set, _after_ removing features not
#'                       present in `ranks` and filtering multimers.}
#' By default, the data frame is sorted by the BH-adjusted p-value and
#' normalized enrichment score (NES).
#' @author Amanda Hiser
#' @references
#'   Subramanian, Tamayo, et al. Gene set enrichment analysis: A knowledge-based
#'   approach for interpreting genome-wide expression profiles.
#'   Proc Natl Acad Sci USA. 102(43):15545-50 (2005).
#'   https://doi.org/10.1073/pnas.0506580102.
#'
#'   Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1alpha-responsive genes
#'   involved in oxidative phosphorylation are coordinately downregulated in
#'   human diabetes. Nat Genet 34, 267–273 (2003).
#'   https://doi.org/10.1038/ng1180.
#' @examples
#' # Prepare ranked input, resolve non-1:1 mapping and split heterodimers
#' ranks <- prepareRanks(stats = t_tests$t_stat,
#'                       features = t_tests$EntrezGeneSymbol)
#'
#' # Run GSEA with defaults
#' go_bp_res <- somaPrGSEA(ranks = ranks)
#' head(go_bp_res$results) # GSEA results
#' head(go_bp_res$final_ranks) # Modified vector used for GSEA
#'
#' # GSEA can be performed using SomaScan AptNames, instead of genes
#' ranks_apt <- t_tests$t_stat
#' names(ranks_apt) <- t_tests$AptName
#' meta <- SomaDataIO::getAnalyteInfo(example_data_11k)
#'
#' # Must set 'use_aptnames = TRUE'
#' res_apt <- somaPrGSEA(ranks = ranks_apt,
#'                       col_meta_df = meta,
#'                       use_aptnames = TRUE)
#' head(res_apt$results)
#'
#' \dontrun{
#' # Using a custom gene set from a GMT file
#' gmt_file <- system.file("extdata", "msigdb_c2_reactome.gmt", package = "SomaEnrich")
#' gmt_list <- fgsea::gmtPathways(gmt_file)
#' gsea_cust <- somaPrGSEA(ranks = rank_vec, cust_paths = gmt_list)
#' }
#' @importFrom fgsea fgseaMultilevel
#' @importFrom SomaDataIO is.apt is.AptName
#' @importFrom withr with_seed
#' @export
somaPrGSEA <- function(ranks,
                       resource = c("bp", "mf", "h", "c1", "c2", "c3",
                                    "c4", "c6", "c7", "c8"),
                       split_ids = TRUE,
                       resolve_multimapping = TRUE,
                       resolve_method = c("abs", "min", "max", "rank"),
                       id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
                       use_aptnames = FALSE,
                       col_meta_df = NULL,
                       cust_paths = NULL,
                       min_feats = 5L,
                       max_feats = 500L,
                       seed = 42,
                       n_perm = 1000L,
                       verbose = interactive(), ...) {

    resource <- match.arg(resource)
    resource <- ifelse(!is.null(cust_paths), "custom", resource)
    resolve_method <- match.arg(resolve_method)

    id_type  <- match.arg(id_type)
    id_col <- switch(id_type, 
                     EntrezGeneID = "entrez_id", 
                     EntrezGeneSymbol = "gene_symbol")

    stopifnot(
        "`ranks` must be a named vector." = !is.null(names(ranks)),
        "`ranks` appears to be SomaScan identifiers, did you set `use_aptnames = TRUE`?" = !(any(SomaDataIO::is.apt(names(ranks))) & !use_aptnames & is.null(cust_paths)),
        "`ranks` must be comprised of a numeric ranking statistic." = is.numeric(unname(ranks))
    )

    # Pre-process ranks by resolving non-1:1 mapping and breaking up heterodimers
    if ( split_ids || resolve_multimapping ) {
        final_ranks <- prepareRanks(stats = unname(ranks), 
                                    features = names(ranks),
                                    split_ids = split_ids,
                                    resolve_multimapping = resolve_multimapping,
                                    resolve_method = resolve_method,
                                    verbose = verbose)
    } else {
        final_ranks <- ranks
    }

    # Prepare pathways for input into fgsea
    if ( !is.null(cust_paths) ) {
        stopifnot("`cust_paths` must be a list of named vectors." = inherits(cust_paths, "list"))
        if ( use_aptnames ) {
            warning("`use_aptnames` is ignored when `cust_paths` is provided. ",
                    "Custom pathways must already use the same identifiers as `ranks`.",
                    call. = FALSE)
        }
        sel_paths <- cust_paths
    } else {
        # Filter pathway map to only retain selected resource
        sel_paths <- pathway_map[pathway_map$group_code == resource, ]
        
        # Add pathway name from addl pathway info lookup table
        name_map <- setNames(path_name_lookup$pathway_name, path_name_lookup$pathway_id)
        sel_paths$pathway_name <- unname(name_map[sel_paths$pathway_id])

        # Convert to list of named vectors for input into fgsea
        sel_paths <- df2list(sel_paths, "pathway_name", id_col, clean_names = FALSE)

        if ( use_aptnames ) {
            stopifnot(
                "AptNames must be provided as input when `use_aptnames = TRUE`." = SomaDataIO::is.AptName(names(ranks)),
                "`col_meta_df` must be provided when `use_aptnames = TRUE`." = !is.null(col_meta_df)
            )

            # Convert gene paths to AptNames.
            # Gene ID type doesn't matter, and is hardcoded as symbols here
            sel_paths <- genePath2aptPath(sel_paths,
                                          col_meta_df = col_meta_df,
                                          id_type = "EntrezGeneSymbol",
                                          verbose = verbose)
        }
    }

    # Stops if names of ranks vector aren't found in chosen/provided pathways,
    # possibly due to ID type mismatch (AptNames for ranks and genes in pathways)
    if ( !any(names(final_ranks) %in% unique(unlist(sel_paths))) ) {
        stop("No elements in `ranks` found in the selected gene sets.")
    }

    path_sizes <- data.frame(pathway = names(sel_paths),
                             starting_set_size = lengths(sel_paths))

    # Determine value to use for fgseaMultilevel(scoreType),
    # see https://github.com/alserglab/fgsea/issues/87
    if ( all(final_ranks < 0) ) {
        score_type <- "neg"
    } else if ( all(final_ranks > 0) ) {
        score_type <- "pos"  # When abs value ranks are provided
    } else {
        score_type <- "std"
    }

    if ( verbose ) {
        writeLines(.rule("Performing Enrichment Analysis", line_col = "cyan"))
        .todo("Running pre-ranked GSEA...")
        if ( split_ids ) {
            .info("If {.code split_ids = TRUE}, any heterodimers present were split into individual elements.")
            .info("This may introduce ties. See {.code ?somaPrGSEA()}.")
        }
    }

    # TODO: Consider `nproc = 1` for reproducibility, per fgsea developers
    results <- withr::with_seed(seed,
                                fgsea::fgseaMultilevel(pathways = sel_paths,
                                                       stats    = final_ranks,
                                                       minSize  = min_feats,
                                                       maxSize  = max_feats,
                                                       nPermSimple = n_perm,
                                                       scoreType = score_type,
                                                       ...)
    )

    # Rearrange/rename cols & merge in path size info
    colnames(results)[colnames(results) == "size"] <- "final_set_size"
    results <- merge(results, path_sizes, by = "pathway")
    col_order <- c("pathway", "pval", "padj", "log2err", "ES", "NES",
                   "leadingEdge", "starting_set_size", "final_set_size")
    results <- results[, col_order]

    # Order by normalized enrichment score & adj p-val
    results <- results[order(results$padj, -abs(results$NES), decreasing = FALSE), ]

    # Add pathway ID
    id_map  <- setNames(path_name_lookup$pathway_id, path_name_lookup$pathway_name)
    new_col <- data.frame(resource_code = resource,
                          pathway_id    = unname(id_map[results$pathway]))

    results <- cbind(new_col, results)
    
    if ( verbose ) {
        .done("Done!")
    }
    
    return(list(results = results,
                final_ranks = final_ranks))
}

