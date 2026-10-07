# Perform Overrepresentation Analysis for SomaScan

Implements the algorithm from
[`fgsea::fora()`](https://rdrr.io/pkg/fgsea/man/fora.html) to perform
overrepresentation analysis (ORA) via hypergeometric test, tailored
specifically to proteomic SomaScan data.

## Usage

``` r
somaORA(
  features,
  universe,
  resource = c("bp", "mf", "h", "c1", "c3", "c4", "c6", "c7", "c8"),
  split_ids = TRUE,
  unique_features = TRUE,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  use_aptnames = FALSE,
  col_meta_df = NULL,
  cust_paths = NULL,
  min_feats = 5L,
  max_feats = 500L,
  verbose = interactive()
)
```

## Arguments

- features:

  Character. A vector containing features of interest, often
  differentially expressed genes or SomaScan analytes. Vector elements
  must be either Entrez gene symbols/identifiers (unless custom gene
  sets are used, see `cust_paths`) or SomaScan identifiers in
  R-compatible `AptName` format (e.g. `seq.1234.56`). If using
  `AptNames`, you must set `use_aptnames = TRUE`.

- universe:

  Character. A vector of unique gene symbols or `AptNames`. This vector
  represents the overall "universe" from which the elements in
  `features` were identified. The `universe` should encompass *all*
  elements in the differential expression analysis from which the set of
  over-represented targets (`features`) was chosen. See `Details` for
  more information. If using `AptNames`, you must set
  `use_aptnames = TRUE`.

- resource:

  Character. String specifying the gene set database resource to use.
  Choose one of:

  "bp"

  :   GO Biological Process

  "mf"

  :   GO Molecular Function

  "h"

  :   MSigDB Hallmark

  "c1"

  :   MSigDB Positional

  "c2"

  :   MSigDB Curated (note: does not contain KEGG or BioCarta)

  "c3"

  :   MSigDB Regulatory Target

  "c4"

  :   MSigDB Computational (note: does not contain CM submodule)

  "c6"

  :   MSigDB Oncogenic Signature

  "c7"

  :   MSigDB Immunologic Signature

  "c8"

  :   MSigDB Cell Type Signature

  The default is `bp`. For more advanced filtering, or to use gene sets
  not provided here, see the `cust_paths` argument. Note that `resource`
  is ignored when a data set is supplied to `cust_paths`.

- split_ids:

  Logical. If the identifiers used in the `features` and `universe`
  vectors contain heterodimers (elements separated by a delimiter, like
  `"GENE1|GENE2"`), should the heterodimers be split into individual
  components, like `c("GENE1", "GENE2")`? If TRUE (the default), this
  may introduce duplicate identifiers into `features`/`universe`.
  Applies to *both* `features` and `universe`. See `Details` and
  examples.

- unique_features:

  Logical. If present, should duplicates in `features` be removed?
  Default is FALSE. If TRUE, the resultant vector will only contain
  unique elements. The duplicated features are often the result of
  setting `split_ids = TRUE` when heterodimers are present in the
  `features` vector. See examples.

- id_type:

  Character. Type of gene identifier used in `features`. Options are
  "EntrezGeneSymbol" or "EntrezGeneID". Default is "EntrezGeneSymbol".
  Ignored if `use_aptnames = TRUE`.

- use_aptnames:

  Logical. Should analysis be performed in aptamer-space, i.e. using
  SomaScan `AptNames` instead of genes? Default is FALSE. If TRUE, the
  input ranking data must contain `AptNames` instead of gene
  identifiers. If using `cust_paths` when `use_aptnames = TRUE`, custom
  gene set members must also be `AptNames`.

- col_meta_df:

  Optional. Data frame of column metadata, likely created by calling
  [`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html)
  on the `soma_adat` used for differential expression analysis. Required
  when `use_aptnames = TRUE`. This object will be used to map between
  SomaScan `AptNames` and gene symbols.

- cust_paths:

  Optional. A named list containing character vectors of gene sets or
  functional groups of interest. If provided, `cust_paths` will override
  `resource`. The name of each list element should be the name of the
  pathway/gene set. The identifier used in each vector should be the
  same type as the names of `ranks`. This object is still subject to the
  minimum/maximum feature cutoffs defined by `min_feats` and
  `max_feats`.

  This argument is meant to be used when the desired biological network
  is not present in the `resource` options.

- min_feats:

  Numeric. Minimum number of features required for a gene set to be
  retained. Default is 5. Sets smaller than this value will be
  discarded. If `NULL`, no minimum is used. See `Details` for more
  information.

- max_feats:

  Numeric. Maximum number of features required for a gene set to be
  retained. Default is 500. See `Details` for more information.

- verbose:

  Logical. Should progress messages be printed to the console? By
  default, messages will only be shown in interactive R sessions. Set to
  TRUE to show all messages, and FALSE to silence messages.

## Value

A `data.frame` object containing the results of ORA, where the results
in each row correspond to a single tested pathway. The data frame
contains the following columns:

- resource_code:

  Abbreviated character string representing the name of the original
  gene set resource.

- pathway_id:

  Pathway identifier/accession number from the original gene set or
  pathway resource.

- pathway:

  Name of the pathway.

- pval:

  An enrichment p-value from hypergeometric test.

- padj:

  The BH-adjusted p-value.

- foldEnrichment:

  The degree of enrichment, relative to the background.

- overlap:

  The size of the overlap between the pathway and features of interest.

- size:

  The size of the pathway.

- overlapFeatures:

  A vector containing the overlapping features (see 'overlap' column).

By default, the data frame is sorted by the BH-adjusted p-value and fold
enrichment.

## Additional Details

Overrepresentation analysis is a statistical method that determines
whether sets of features, often genes, are present in a data set more
than would be expected by chance. As an example, imagine a SomaScan
experiment in which 200 analytes are determined to be differentially
expressed in a group of breast cancer samples. Using ORA, a researcher
can investigate whether any pathways, like cellular proliferation, are
more highly represented in those 200 analytes than would be expected by
chance.

## Feature Selection

The primary input to ORA is a set of biological features of interest.
Selecting these features typically involves defining a specific subset,
often differentially expressed genes, based on previously defined
thresholds or cutoffs. For example, one could use the combined
statistical significance and log fold change (logFC) cutoffs of p \<
0.05 and logFC \> 2, respectively, to identify significantly
differentially expressed genes. This list would then be used as input
for ORA. Please note that the most appropriate cutoff(s) to use will
likely vary by experiment; this is only an example.

## Enrichment Universe

The `universe` (sometimes referred to as "background") serves as a point
of reference against which the genes or analytes of interest (provided
to `features`) can be compared. Setting the correct background universe
is *crucial* for accurate ORA results. The universe must encompass all
features in the analysis results (typically differential expression
analysis). If analyte filtering was performed prior to differential
expression analysis, the filtered analytes should not be used to set the
ORA universe.

If `features` contains elements that are not in `universe`, those
elements will be dropped from `features`. Additionally, `universe`
elements that are not found in any gene set will be dropped.

## Handling Heterodimers

When a SOMAmer Reagent targets a protein heterodimer, the genes
associated with that heterodimer are collapsed using a delimiter
character (ex. "CGA\|CGB3\|CGB7") for concise display in the SomaScan
menu. However, when performing ORA, these collapsed and delimited
strings must be split into individual gene symbols, so each symbol can
be matched against gene-based databases like GO or MSigDB. In
SomaEnrich, the `resource` or `cust_paths` arguments direct annotation
of the features of interest to known pathways or gene sets. If the
elements of `features` do not match elements in `resource`/`cust_paths`,
the non-matching elements in `features`/`universe` will be dropped.
Consequently, to preserve these elements, the default behavior is to
split delimited strings.

## Min/Max Set Sizes

The statistical power of GSEA/ORA to detect enrichment depends on the
size of the gene set being tested, and GSEA normalization is not
accurate for very small feature sets. The minimum allowable gene set is
defined by the `min_feats` argument; any sets smaller or larger than
`min_feats` or `max_feats` will be discarded *before* performing
GSEA/ORA. Please note that this threshold is applied *after* filtering
to remove features from each gene set that are not also present in the
`ranks`, i.e. taking the intersect of the `ranks` and gene lists. The
default `min_feats` threshold used here differs from other enrichment
analysis tools. The 15-gene minimum suggested by the Broad Institute's
[GSEA software](https://www.gsea-msigdb.org/gsea/index.jsp) was
established in a whole-genome context, where the previously described
intersection step removes few genes. However, in the reduced feature
space of SomaScan, more genes are lost from taking this intersection. In
some cases, this can greatly reduce the size of available gene sets. To
account for this loss, a default minimum threshold of 5 is used.
Consider final gene set size accordingly when reviewing your enrichment
analysis results.

## Author

Amanda Hiser

## Examples

``` r
# Genes mapping to differentially expressed analytes
deg <- t_tests[t_tests$fc > 1 & t_tests$fdr < 0.05, ]$EntrezGeneSymbol

# Input universe = all genes in statistical test results
uni <- t_tests$EntrezGeneSymbol

# To properly check for enrichment, values in the input vector must match 
# values in the gene sets/pathway resource. 
# SomaScan contains heterodimers that won't match non SomaScan-derived gene sets
head(deg)
#> [1] "PZP"           "CGA|FSHB"      "ENPP2"         "CGA|LHB"      
#> [5] "CGA|CGB3|CGB7" "NTM"          
deg[grepl("\\|", deg)]
#> [1] "CGA|FSHB"      "CGA|LHB"       "CGA|CGB3|CGB7" "LTA|LTB"      

# Splitting these values can create duplicates in the input vector
split_deg <- strsplit(deg, "\\|") |> unlist()
split_deg[duplicated(split_deg)]
#> [1] "CGA"     "CGA"     "LEP"     "NTM"     "LHB"     "SHBG"   
#> [7] "MGAT5"   "C1QTNF4" "ROBO2"  

# The following will run ORA with GO BP gene sets (the default),
# split heterodimers into individual components (the default),
# and remove duplicates created from splitting
res <- somaORA(features = deg,
               universe = uni,
               unique_features = TRUE)

# ORA can be performed using SomaScan AptNames, instead of genes
dea <- head(t_tests$AptName, 50)
apt_uni <- t_tests$AptName

# A progress bar will appear to track AptName conversion
meta <- SomaDataIO::getAnalyteInfo(example_data_11k)
res_apt <- somaORA(features = dea,
                   use_aptnames = TRUE,
                   universe = apt_uni,
                   col_meta_df = meta)
```
