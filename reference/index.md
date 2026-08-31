# Package index

## Enrichment Analysis

Core functions for performing pathway enrichment analyses with SomaScan
data.

- [`somaORA()`](https://somalogic.github.io/SomaEnrich/reference/somaORA.md)
  : Perform Overrepresentation Analysis for SomaScan
- [`somaPrGSEA()`](https://somalogic.github.io/SomaEnrich/reference/somaPrGSEA.md)
  : Perform Pre-Ranked GSEA for SomaScan
- [`prepareRanks()`](https://somalogic.github.io/SomaEnrich/reference/prepareRanks.md)
  : Prepare Ranking Metrics for SomaScan-Based Enrichment Analyses

## ID Mapping

Functions for converting between SomaScan analyte identifiers and gene
identifiers.

- [`apt2gene()`](https://somalogic.github.io/SomaEnrich/reference/apt2gene.md)
  : Convert SomaScan Analytes to Genes
- [`gene2apt()`](https://somalogic.github.io/SomaEnrich/reference/gene2apt.md)
  : Convert Gene Identifiers to SomaScan Analytes (v2)
- [`go2apt()`](https://somalogic.github.io/SomaEnrich/reference/go2apt.md)
  : Retrieve SomaScan Analytes in a Given GO Term
- [`collapseAptData()`](https://somalogic.github.io/SomaEnrich/reference/collapseAptData.md)
  : Transform Aptamer-Centric Data to Gene-Centric Data

## Pathway Conversion

Functions for converting pathway/feature set representations between
gene-based and SomaScan analyte-based formats.

- [`aptPath2genePath()`](https://somalogic.github.io/SomaEnrich/reference/aptPath2genePath.md)
  : Convert SomaScan ID-Based Feature Set into Gene Set
- [`genePath2aptPath()`](https://somalogic.github.io/SomaEnrich/reference/genePath2aptPath.md)
  : Convert Gene Set into SomaScan ID-Based Feature Set
- [`gmt_gene2apt()`](https://somalogic.github.io/SomaEnrich/reference/gmt_gene2apt.md)
  : Convert Gene-Based GMT File to AptNames

## Plotting

Visualization functions for enrichment analysis results.

- [`plotBubble()`](https://somalogic.github.io/SomaEnrich/reference/plotBubble.md)
  : ORA Bubble Plot
- [`plotES()`](https://somalogic.github.io/SomaEnrich/reference/plotES.md)
  : GSEA Enrichment Plot with Leading Edge Annotation

## Utilities

Helper functions for data preparation and manipulation.

- [`df2list()`](https://somalogic.github.io/SomaEnrich/reference/df2list.md)
  : Convert a Data Frame to a Named List of Vectors

## Data Objects

Objects provided with `SomaEnrich`.

- [`example_data_11k`](https://somalogic.github.io/SomaEnrich/reference/example_data_11k.md)
  : Example 11K Data Set
- [`t_tests`](https://somalogic.github.io/SomaEnrich/reference/t_tests.md)
  : t-test Example Data
- [`pathway_map`](https://somalogic.github.io/SomaEnrich/reference/pathway_map.md)
  : Biological Networks for Enrichment Analysis
- [`ex_apt_pathway`](https://somalogic.github.io/SomaEnrich/reference/ex_apt_pathway.md)
  : Example Aptamer-Based Pathway
- [`ex_gene_pathway`](https://somalogic.github.io/SomaEnrich/reference/ex_gene_pathway.md)
  : Example Gene-Based Pathway
