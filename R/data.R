#' t-test Example Data
#'
#' Results of performing differential expression (via t-test) with an example 
#' SomaScan dataset (`example_data_11k`). The `Sex` variable was used to compare 
#' males (N = 81) against females (N = 81). The results are sorted by 
#' t-statistic and adjusted p-value.
#' 
#' The SomaDataIO "Two-Group Comparison" workflow was used to produce this data 
#' set, with modifications (different input ADAT, and fold change and 
#' log2 fold change values were added to the final table).
#' @name t_tests
#' @docType data
#' @keywords datasets
#' @source \url{https://somalogic.github.io/SomaDataIO/articles/stat-two-group-comparison.html}
#' @examples
#' colnames(t_tests)
#'
#' head(t_tests)
#'
#' # Formula used to generate results for a given analyte
#' t_tests$formula[[1]]
NULL


#' Biological Networks for Enrichment Analysis
#'
#' @description
#' The `pathway_map` data object is a `data.frame` that maps
#' genes to gene sets/pathways from public data repositories.
#'
#' In the `pathway_map` object, the columns are:
#' \itemize{
#'   \item `gene_symbol`: Gene symbol (character)
#'   \item `entrez_id`: Entrez gene identifier (character)
#'   \item `pathway_id`: Pathway identifier from respective data source (character)
#'   \item `group_code`: Abbreviated code for the data
#'     source/collection combination, e.g. `"bp"`, `"h"`, or `"c2"` (factor)
#' }
#'
#' The following networks are available, identified by their `group_code`:
#' \itemize{
#'   \item `bp` (GO Biological Process): describes the larger cellular or
#'   physiological role played by a gene product, often in coordination with
#'   other genes
#'   \item `mf` (GO Molecular Function): describes activities performed by
#'   gene products at the molecular level, such as binding or catalysis
#'   \item `h` (MSigDB Hallmark): summarizes well-defined biological states
#'   \item `c1` (MSigDB Positional): gene sets corresponding to human
#'   chromosome cytogenetic bands
#'   \item `c2` (MSigDB Curated): selected gene sets from online pathway
#'   databases and biomedical literature (KEGG and BioCarta removed for
#'   licensing reasons)
#'   \item `c3` (MSigDB Regulatory Target): potential targets of regulation
#'   by transcription factors or miRNAs
#'   \item `c4` (MSigDB Computational): defined by mining cancer-oriented
#'   expression data (KEGG pathways removed for licensing reasons)
#'   \item `c6` (MSigDB Oncogenic Signature): signatures of cellular pathways
#'   often dysregulated in cancer
#'   \item `c7` (MSigDB Immunologic Signature): gene sets associated with
#'   immune disease
#'   \item `c8` (MSigDB Cell Type Signature): curated markers for cell types
#'   identified in single-cell sequencing studies of human tissue
#' }
#' 
#' @details
#' [Gene Ontology](https://geneontology.org/) data (retrieved via the `GO.db`
#' R package, v3.23.1) is made available under the terms of the 
#' [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).
#' 
#' [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) data (version 2026.1.Hs,
#' retrieved via the `msigdbr` R package v26.1.0) is copyright
#' © 2004–2025 Broad Institute, Inc., Massachusetts Institute of Technology,
#' and Regents of the University of California, and is made available under
#' the terms of the 
#' [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).
#' 
#' Data was retrieved from these resources on 06/01/2026.
#' @name pathway_map
#' @docType data
#' @author Amanda Hiser
#' @keywords datasets
#' @examples
#' head(pathway_map)
#' length(unique(pathway_map$pathway_id))
#' table(pathway_map$group_code)
NULL


#' Example Gene-Based Pathway
#'
#' Example biological pathway comprised of gene identifiers (Entrez Gene IDs). 
#' Pathway is a subset of "HALLMARK_INFLAMMATORY_RESPONSE" (M5932) from MSigDB.
#' @name ex_gene_pathway
#' @docType data
#' @keywords datasets
#' @source \url{gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INFLAMMATORY_RESPONSE}
#' @examples
#' head(ex_gene_pathway, n = 20)
#' 
#' length(ex_gene_pathway)
NULL


#' Example Aptamer-Based Pathway
#'
#' Example biological pathway comprised of SomaScan analyte identifiers 
#' (`AptNames`). Pathway is a subset of "HALLMARK_INFLAMMATORY_RESPONSE" 
#' (M5932) from MSigDB, converted to `AptNames` via `gene2apt()`. Only genes
#' mapping to an analyte in the SomaScan menu were retained.
#' @name ex_apt_pathway
#' @docType data
#' @keywords datasets
#' @source \url{gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INFLAMMATORY_RESPONSE}
#' @examples
#' head(ex_apt_pathway, n = 20)
#' 
#' length(ex_apt_pathway)
NULL


#' Example 11K Data Set
#' 
#' The `example_data_11k` object is intended to provide existing and prospective
#' SomaLogic customers with example data to enable analysis preparation prior
#' to receipt of SomaScan data, and also for those generally curious about the
#' SomaScan data deliverable. It is **not** intended to be used as a control
#' group for studies or provide any metrics for SomaScan data in general.
#' @name example_data_11k
#' @docType data
#' @keywords datasets
#' @source \url{https://github.com/SomaLogic/SomaLogic-Data/blob/main/example_data_v5.0_plasma.adat}
#' @examples
#' # S3 print method
#' example_data_11k
#'
#' # Print header info
#' print(example_data_11k, show_header = TRUE)
#'
#' # View object class
#' class(example_data_11k)
#'
#' # Retrieve annotations
#' SomaDataIO::getAnalyteInfo(example_data_11k)
NULL
