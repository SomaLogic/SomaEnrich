# Perform Pre-Ranked GSEA for SomaScan

Wrapper function for
[`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) to perform
pre-ranked gene set enrichment analysis (GSEA), tailored specifically to
proteomic SomaScan data.

## Usage

``` r
somaPrGSEA(
  ranks,
  resource = c("bp", "mf", "h", "c1", "c2", "c3", "c4", "c6", "c7", "c8"),
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
  verbose = interactive(),
  ...
)
```

## Arguments

- ranks:

  Either the output of
  [`prepareRanks()`](https://somalogic.github.io/SomaEnrich/reference/prepareRanks.md),
  or a named vector of a numeric ranking statistic. This vector will be
  used to rank features, according to the value of the statistic. Any
  numeric metric may be used (t-statistic, -log10(p), fold change,
  etc.). The vector *must* be named using Entrez gene
  symbols/identifiers *or* SomaScan identifiers in R-compatible
  `AptName` format (e.g. `seq.1234.56`). If using `AptNames`, set
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

  Logical. If the identifiers used for the statistical metric contain
  heterodimers (features separated by a delimiter, like
  `"GENE1|GENE2"`), should the heterodimers be split into individual
  components, like `c("GENE1", "GENE2")`? If TRUE (the default), the
  final ranking metric will contain duplicated numeric values, one for
  each unique member of the heterodimer. See `Details` and examples.

- resolve_multimapping:

  Logical. Should many-to-one relationships between SOMAmer Reagents and
  genes (i.e. multiple SOMAmer Reagents mapping to the same gene) in the
  statistical ranking metric be resolved? Default is TRUE. Only one
  entry per gene will be retained, based on the `resolve_method`
  selection criteria. See `Details`.

- resolve_method:

  Character. Selection method for choosing a single representative value
  for a gene when non-1:1 mapping exists between analytes and genes.
  Choose one of:

  "abs"

  :   (Default) Selects the entry with the largest absolute value per
      gene. Recommended for signed statistics (e.g. t-statistic, log
      fold change) in bidirectional analyses — unlike `"rank"` or
      `"max"`, it preserves strong signal symmetrically at both the top
      and bottom of the ranking vector.

  "min"

  :   Selects the entry with the minimum numeric value, including
      negative values. Use for metrics where smaller values are better
      (e.g., p-values) or when smaller/negative signed values are
      important.

  "max"

  :   Selects the entry with the maximum numeric value. Use for metrics
      where larger values are preferred (e.g. positive one-sided tests).

  "rank"

  :   Selects the entry with the smallest positional index (i.e. first
      occurrence in the input vector). Use when input is pre-sorted by
      one or more variables, and the "best" entry is the one closest to
      position 1. This means that the selected entry depends on the
      order of the input, so `ranks` must be **pre-sorted** to place
      values of highest importance first.

- id_type:

  Character. Type of gene identifier used in `ranks`. Options are
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

- seed:

  Numeric. Seed to set for GSEA calculations. Because the GSEA algorithm
  used here is non-deterministic, a seed must be set for results to be
  reproducible. Default is 42.

- n_perm:

  Numeric. Number of permutations to perform. Default is 1000.

- verbose:

  Logical. Should progress messages be printed to the console? By
  default, messages will only be shown in interactive R sessions. Set to
  TRUE to show all messages, and FALSE to silence messages.

- ...:

  Optional arguments passed to
  [`fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html).

## Value

A list with two elements:

- results:

  A `data.frame` object containing the results of GSEA, where the
  results in each row correspond to a single tested gene set.

- final_ranks:

  The final ranks used as input for GSEA. If
  `resolve_multimapping = TRUE` and `split_ids = TRUE`, the ranks used
  as input were modified to split protein heterodimers into unique gene
  members and resolve non-1:1 mapping. Therefore, this vector may differ
  from what was provided as input to `somaPrGSEA()`.

The `results` data frame contains the following columns:

- resource_code:

  Abbreviated character string representing the name of the original
  gene set resource.

- pathway_id:

  Pathway identifier/accession number from the original gene set or
  pathway resource.

- pathway:

  Name of the pathway.

- pval:

  Enrichment p-value.

- padj:

  FDR-adjusted p-value.

- log2err:

  The expected error for the standard deviation of the p-value
  logarithm.

- ES:

  The enrichment score calculated via the [Broad GSEA
  implementation](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Enrichment_Score_(ES)).
  The ES represents the degree to which a set `S` is over-represented at
  the top or bottom of a ranked list of genes (`L`). The ES is
  calculated by walking down the list `L`, increasing a running-sum
  statistic whenever a gene in `L` is found in set `S`, and decreasing
  it each time a gene is not found in `S`. The final ES is the maximum
  deviation from zero encountered in the random walk. A positive ES
  indicates that `S` is enriched toward the top of the ranked list,
  while a negative ES indicates that `S` is enriched at the bottom of
  the ranked list.

- NES:

  Enrichment score normalized to account for the size of each gene set,
  calculated by taking the ES and dividing it by the mean ES for all
  data set permutations. The NES is the primary statistic to use for
  examining GSEA results.

- leadingEdge:

  Vector with indexes of leading edge features that drive the
  enrichment, see [Running a Leading Edge
  Analysis](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading).
  Leading edge features in the set `S` appear in the ranked list `L` at,
  or before, the point where the running sum reaches its maximum
  deviation from zero. Essentially, the leading edge subset is the core
  of a feature set that accounts for the enrichment signal.

- starting_set_size:

  Original size of the pathway/gene set, without modifications.

- final_set_size:

  Number of SomaScan analytes mapping to features in the set, *after*
  removing features not present in `ranks` and filtering multimers.

By default, the data frame is sorted by the BH-adjusted p-value and
normalized enrichment score (NES).

## Additional Details

Gene Set Enrichment Analysis (GSEA) is a rank-based statistical method
that tests whether predefined groups of biologically related entities
are over-represented at the top or bottom of a ranked list of features.
Unlike many other methods of biological interpretation (like ORA), GSEA
does not rely on arbitrary significance cutoffs. In proteomics
experiments, GSEA is used to detect coordinated changes across groups of
functionally related proteins/genes, and provides a source for
functional biological interpretation of statistical test results.

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

## SOMAmer Reagent and Target Relationships

Non-1:1 mapping between SOMAmer Reagents and their targets is a known
characteristic of the SomaScan menu. This can lead to a single gene
being associated with multiple statistical metric results. For example,
prostate-specific antigen (PSA) and its associated gene, KLK3, are
targeted by 4 SOMAmer Reagents in the 7K SomaScan menu. This means that,
*at the gene level*, KLK3 has 4 RFU values in any 7K SomaScan data set.

In pre-ranked GSEA, however, values in the input ranking data must be
unique. Therefore, the input data cannot contain 4 measurements for
PSA - only one must be chosen to be representative of the KLK3 gene.
There are many acceptable selection methods, and the most appropriate
one will depend on your question and/or experimental design. See the
analyte-to-gene mapping vignette
(`browseVignettes(package = "SomaEnrich")`) for more details.

## References

Subramanian, Tamayo, et al. Gene set enrichment analysis: A
knowledge-based approach for interpreting genome-wide expression
profiles. Proc Natl Acad Sci USA. 102(43):15545-50 (2005).
https://doi.org/10.1073/pnas.0506580102.

Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1alpha-responsive
genes involved in oxidative phosphorylation are coordinately
downregulated in human diabetes. Nat Genet 34, 267–273 (2003).
https://doi.org/10.1038/ng1180.

## Author

Amanda Hiser

## Examples

``` r
# Prepare ranked input, resolve non-1:1 mapping and split heterodimers
ranks <- prepareRanks(stats = t_tests$t_stat,
                      features = t_tests$EntrezGeneSymbol)

# Run GSEA with defaults
go_bp_res <- somaPrGSEA(ranks = ranks)
#> Warning: There are ties in the preranked stats (0.44% of the list).
#> The order of those tied genes will be arbitrary, which may produce unexpected results.
head(go_bp_res$results) # GSEA results
#>   resource_code pathway_id                    pathway         pval
#> 1            bp GO:0007416           synapse assembly 2.706503e-07
#> 2            bp GO:0050808       synapse organization 3.435995e-07
#> 3            bp GO:0034330 cell junction organization 1.746782e-07
#> 4            bp GO:0008038         neuron recognition 4.620156e-07
#> 5            bp GO:0099054        presynapse assembly 5.834195e-07
#> 6            bp GO:0034329     cell junction assembly 9.945292e-07
#>           padj   log2err        ES      NES
#> 1 0.0006940710 0.6749629 0.4670968 1.997156
#> 2 0.0006940710 0.6749629 0.3987822 1.838772
#> 3 0.0006940710 0.6901325 0.3669554 1.728795
#> 4 0.0006999536 0.6749629 0.7133333 2.342164
#> 5 0.0007071045 0.6594444 0.6672418 2.238701
#> 6 0.0010044745 0.6435518 0.3986414 1.799007
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     leadingEdge
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                NTM, ROBO2, SLITRK1, LSAMP, RELN, GHRL, OPCML, NRCAM, NLGN1, C1QL3, CLSTN1, MDGA1, PTPRD, SLITRK5, C1QL2, SLITRK6, LRRTM2, EPHB2, CRTAC1, GRID2, IGLON5, NPTX1, CDH1, LRRC15, NRXN3, NLGN4Y, PLXNA1, EFNB2, CBLN1, THBS2, HAPLN4, LRRTM4, IL1RAP, SLITRK2, CSMD2, EGFLAM, NLGN2, NEGR1, IL1RAPL2, LRFN4, ABL1, UBE2M, FGFR1, LRRC4, LRFN3, DAG1, NRXN1, NLGN3, WNT7A, LRRTM3, LRP4, NECTIN1, ICAM5, ADGRL3, KIRREL3, LRFN1, RAP2A, NRG3, LIN7B
#> 2                                                                                                                                                       NTM, ROBO2, SLITRK1, LSAMP, RELN, SEZ6L, GHRL, OPCML, GRN, NRCAM, NLGN1, C1QL3, CLSTN1, MDGA1, PTPRD, SLITRK5, C1QL2, SLITRK6, CHRDL1, LRRTM2, EPHB2, CAST, ITGB3, CRTAC1, GRID2, LRP8, IGLON5, SPARCL1, NPTX1, CDH1, LRRC15, ERBB3, NRXN3, CNTN4, NLGN4Y, PLXNA1, EFNB2, CBLN2, CBLN1, NTNG1, TREM2, FCGR2B, APBB2, DIXDC1, C1QL1, C3, THBS2, HAPLN4, CNTN6, LRRTM4, NRP2, IL1RAP, SLITRK2, CSMD2, EGFLAM, CAMK1, NLGN2, NEGR1, IL1RAPL2, INA, PALM, LRFN4, LAMB1, LAMC1, ABL1, UBE2M, SORT1, FGFR1, PGRMC1, LRRC4, LRFN3, DAG1, TENM3, NRXN1, PTK7, LARGE1, NLGN3, ARHGAP22, WNT7A, CNTN2, LRRTM3, AMOT, LRP4, NECTIN1, ICAM5, ADGRL3, CPNE6, KIRREL3, PRNP, LRFN1, CX3CL1, HIP1R, RAP2A, NRG3, LIN7B
#> 3 NTM, ROBO2, SLITRK1, CNTNAP2, LSAMP, RELN, THY1, SEZ6L, GHRL, OPCML, GRN, DSG2, NRCAM, NLGN1, C1QL3, CLSTN1, MDGA1, PTPRD, AGT, SLITRK5, C1QL2, SLITRK6, CHRDL1, CSF1R, LRRTM2, EPHB2, CAST, ITGB3, CRTAC1, GRID2, LRP8, IGLON5, SPARCL1, GDF2, NPTX1, CDH1, LRRC15, ERBB3, NRXN3, CNTN4, NLGN4Y, PLXNA1, HRG, EFNB2, CBLN2, CBLN1, NTNG1, TREM2, GREM1, FCGR2B, APBB2, DIXDC1, C1QL1, C3, THBS2, HAPLN4, CNTN6, LRRTM4, NRP2, IL1RAP, SLITRK2, LDB1, MICALL2, CSMD2, IL17A, EGFLAM, VCL, CDH4, CAMK1, NLGN2, NEGR1, IL1RAPL2, INA, PALM, LRFN4, LAMB1, LAMC1, MMP14, ANK2, ABL1, MPDZ, UBE2M, SORT1, FGFR1, PGRMC1, LRRC4, LRFN3, CD9, DAG1, TGFB2, TENM3, NRXN1, PTK7, LARGE1, NLGN3, ARHGAP22, WNT7A, KIRREL1, CNTN2, LRRTM3, AMOT, LRP4, NECTIN1, ICAM5, ADGRL3, ARL2, CDH15, CPNE6, KIRREL3, PRNP, LRFN1, CX3CL1, HIP1R, EPHA3, PKN2, RAP2A, NRG3, LIN7B
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 NTM, ROBO2, CNTNAP2, OPCML, NRCAM, EPHB2, CRTAC1, IGLON5, CNTN4, CNTN6, ROBO1
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  SLITRK1, NLGN1, MDGA1, PTPRD, GRID2, NLGN4Y, EFNB2, CBLN1, IL1RAP, SLITRK2, NLGN2, IL1RAPL2, LRFN4, LRFN3, NRXN1, NLGN3, WNT7A, LRRTM3, LRP4
#> 6                                                                                                                                                                                                                                                                                                                                                                                                   NTM, ROBO2, SLITRK1, CNTNAP2, LSAMP, RELN, THY1, GHRL, OPCML, DSG2, NRCAM, NLGN1, C1QL3, CLSTN1, MDGA1, PTPRD, AGT, SLITRK5, C1QL2, SLITRK6, LRRTM2, EPHB2, CRTAC1, GRID2, IGLON5, GDF2, NPTX1, CDH1, LRRC15, NRXN3, NLGN4Y, PLXNA1, HRG, EFNB2, CBLN1, GREM1, THBS2, HAPLN4, LRRTM4, IL1RAP, SLITRK2, LDB1, MICALL2, CSMD2, IL17A, EGFLAM, VCL, CDH4, NLGN2, NEGR1, IL1RAPL2, LRFN4, LAMC1, MMP14, ANK2, ABL1, MPDZ, UBE2M, FGFR1, LRRC4, LRFN3, CD9, DAG1
#>   starting_set_size final_set_size
#> 1               264            162
#> 2               562            309
#> 3               803            416
#> 4                42             33
#> 5                46             39
#> 6               490            264
head(go_bp_res$final_ranks) # Modified vector used for GSEA
#>       PZP       CGA      FSHB     ENPP2       LHB      CGB3 
#> 14.285305  9.528913  9.528913  8.862328  8.709291  7.895289 

# GSEA can be performed using SomaScan AptNames, instead of genes
ranks_apt <- t_tests$t_stat
names(ranks_apt) <- t_tests$AptName
meta <- SomaDataIO::getAnalyteInfo(example_data_11k)

# Must set 'use_aptnames = TRUE'
res_apt <- somaPrGSEA(ranks = ranks_apt,
                      col_meta_df = meta,
                      use_aptnames = TRUE)
head(res_apt$results)
#>   resource_code pathway_id                        pathway         pval
#> 1            bp GO:0007416               synapse assembly 2.120392e-10
#> 2            bp GO:0008038             neuron recognition 9.702891e-10
#> 3            bp GO:0019731 antibacterial humoral response 4.898524e-09
#> 4            bp GO:0050808           synapse organization 4.719000e-09
#> 5            bp GO:0051963 regulation of synapse assembly 6.853876e-09
#> 6            bp GO:0034329         cell junction assembly 1.211029e-08
#>           padj   log2err         ES       NES
#> 1 1.414089e-06 0.8266573  0.4866897  2.132269
#> 2 3.235429e-06 0.7881868  0.7347060  2.472118
#> 3 8.167065e-06 0.7614608 -0.6434427 -2.505318
#> 4 8.167065e-06 0.7614608  0.3984248  1.843070
#> 5 9.141700e-06 0.7614608  0.5078850  2.132170
#> 6 1.346058e-05 0.7477397  0.4106525  1.859923
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  leadingEdge
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  seq.8428.102, seq.20550.38, seq.5116.62, seq.10743.13, seq.2999.6, seq.21713.11, seq.22578.17, seq.23581.131, seq.15622.13, seq.5109.24, seq.15620.4, seq.10907.116, seq.21707.15, seq.15521.4, seq.16900.29, seq.9296.15, seq.4568.17, seq.6423.66, seq.16916.19, seq.6904.14, seq.8225.86, seq.5632.6, seq.12758.47, seq.8447.11, seq.6478.2, seq.9256.78, seq.8052.115, seq.14759.149, seq.6557.50, seq.5111.15, seq.8348.4, seq.21710.16, seq.9005.16, seq.8772.5, seq.9313.27, seq.15539.15, seq.16323.8, seq.3339.33, seq.6455.52, seq.6572.10, seq.14048.7, seq.23000.22, seq.2501.51, seq.21385.5, seq.9971.5, seq.12338.27, seq.2630.12, seq.9772.153, seq.7050.5, seq.5082.51, seq.21696.80, seq.3341.33, seq.19111.10, seq.15553.22, seq.5532.53, seq.21690.31, seq.21691.27, seq.8646.61, seq.8369.102, seq.5110.84, seq.17427.26, seq.13109.82, seq.4889.82, seq.4453.83, seq.19558.10, seq.9300.13, seq.14111.15, seq.8245.27, seq.20578.10
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      seq.8428.102, seq.20550.38, seq.5116.62, seq.6965.19, seq.22578.17, seq.15622.13, seq.5109.24, seq.10907.116, seq.8225.86, seq.5632.6, seq.6478.2, seq.8348.4, seq.3298.52, seq.20561.15, seq.5740.17
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               seq.8468.19, seq.21232.39, seq.13699.6, seq.4330.4, seq.4874.3, seq.4982.54, seq.9250.87, seq.22974.25, seq.22403.13, seq.14143.8, seq.4126.22, seq.4135.84, seq.4916.2, seq.15481.45, seq.9384.17, seq.13671.40, seq.5741.55, seq.15576.158
#> 4 seq.8428.102, seq.20550.38, seq.5116.62, seq.10743.13, seq.2999.6, seq.21713.11, seq.19563.3, seq.22578.17, seq.23581.131, seq.15622.13, seq.4992.49, seq.5109.24, seq.15620.4, seq.10907.116, seq.21707.15, seq.15521.4, seq.16900.29, seq.9296.15, seq.4568.17, seq.6423.66, seq.16916.19, seq.3362.61, seq.6904.14, seq.8225.86, seq.3026.5, seq.20187.10, seq.5632.6, seq.12758.47, seq.8447.11, seq.3323.37, seq.6478.2, seq.4467.49, seq.9256.78, seq.8052.115, seq.14759.149, seq.6557.50, seq.2617.56, seq.5111.15, seq.8348.4, seq.3298.52, seq.21710.16, seq.9005.16, seq.8772.5, seq.21887.2, seq.9313.27, seq.5637.81, seq.15539.15, seq.5635.66, seq.3310.62, seq.12753.6, seq.13441.30, seq.16323.8, seq.6404.20, seq.4480.59, seq.3339.33, seq.6455.52, seq.20561.15, seq.6572.10, seq.15387.44, seq.14048.7, seq.23000.22, seq.2501.51, seq.21385.5, seq.9971.5, seq.12338.27, seq.11851.21, seq.2630.12, seq.3592.4, seq.9772.153, seq.13707.27, seq.7050.5, seq.5082.51, seq.11436.6, seq.18306.1, seq.21696.80, seq.18347.15, seq.16300.4, seq.3341.33, seq.19111.10, seq.15553.22, seq.11300.32, seq.5532.53, seq.6590.54, seq.7863.50, seq.21690.31, seq.21691.27, seq.8646.61, seq.8369.102, seq.11107.25, seq.5110.84, seq.9525.1, seq.7935.26, seq.17427.26, seq.13109.82, seq.24946.9, seq.4900.8, seq.4889.82, seq.3296.92, seq.4453.83, seq.25296.3, seq.19558.10, seq.9300.13, seq.14111.15, seq.8245.27, seq.20578.10, seq.23689.52, seq.4557.61, seq.6545.58, seq.7910.41, seq.2827.23, seq.24960.48, seq.9885.41, seq.10981.56, seq.15612.5
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              seq.8428.102, seq.20550.38, seq.5116.62, seq.10743.13, seq.2999.6, seq.22578.17, seq.23581.131, seq.15622.13, seq.15620.4, seq.10907.116, seq.15521.4, seq.16900.29, seq.9296.15, seq.4568.17, seq.16916.19, seq.6904.14, seq.8225.86, seq.5632.6, seq.12758.47, seq.8447.11, seq.6478.2, seq.8052.115, seq.6557.50, seq.8348.4, seq.9313.27, seq.15539.15, seq.3339.33, seq.6572.10, seq.14048.7, seq.23000.22, seq.21385.5, seq.2630.12, seq.9772.153, seq.7050.5, seq.5082.51, seq.21696.80, seq.19111.10, seq.5532.53, seq.21691.27, seq.8646.61, seq.5110.84, seq.17427.26, seq.13109.82, seq.4889.82, seq.4453.83, seq.19558.10, seq.9300.13, seq.14111.15, seq.8245.27
#> 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        seq.8428.102, seq.20550.38, seq.5116.62, seq.10743.13, seq.6965.19, seq.2999.6, seq.21713.11, seq.20564.53, seq.22578.17, seq.23581.131, seq.15622.13, seq.9484.75, seq.5109.24, seq.15620.4, seq.10907.116, seq.21707.15, seq.15521.4, seq.16900.29, seq.9296.15, seq.3484.60, seq.4568.17, seq.6423.66, seq.16916.19, seq.6904.14, seq.8225.86, seq.5632.6, seq.12758.47, seq.8447.11, seq.6478.2, seq.4880.21, seq.9256.78, seq.8052.115, seq.14759.149, seq.6557.50, seq.5111.15, seq.8348.4, seq.21710.16, seq.9005.16, seq.4996.66, seq.8772.5, seq.9313.27, seq.15539.15, seq.18878.15, seq.16323.8, seq.3339.33, seq.6455.52, seq.6572.10, seq.14048.7, seq.23000.22, seq.2501.51, seq.21385.5, seq.25926.29, seq.12891.1, seq.9971.5, seq.9170.24, seq.12338.27, seq.8750.46, seq.20517.1, seq.20589.5, seq.2630.12, seq.9772.153, seq.7050.5, seq.5082.51, seq.21696.80, seq.18347.15, seq.5002.76, seq.7624.19, seq.3341.33, seq.14036.116, seq.19111.10, seq.15553.22, seq.5532.53, seq.21690.31, seq.21691.27, seq.17449.23, seq.8646.61, seq.8369.102
#>   starting_set_size final_set_size
#> 1               213            213
#> 2                42             42
#> 3                67             64
#> 4               408            403
#> 5               150            150
#> 6               342            342

if (FALSE) { # \dontrun{
# Using a custom gene set from a GMT file
gmt_file <- system.file("extdata", "msigdb_c2_reactome.gmt", package = "SomaEnrich")
gmt_list <- fgsea::gmtPathways(gmt_file)
gsea_cust <- somaPrGSEA(ranks = rank_vec, cust_paths = gmt_list)
} # }
```
