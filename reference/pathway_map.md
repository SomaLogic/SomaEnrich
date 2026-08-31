# Biological Networks for Enrichment Analysis

The `pathway_map` data object is a `data.frame` that maps genes to gene
sets/pathways from public data repositories.

In the `pathway_map` object, the columns are:

- `gene_symbol`: Gene symbol (character)

- `entrez_id`: Entrez gene identifier (character)

- `pathway_id`: Pathway identifier from respective data source
  (character)

- `group_code`: Abbreviated code for the data source/collection
  combination, e.g. `"bp"`, `"h"`, or `"c2"` (factor)

The following networks are available, identified by their `group_code`:

- `bp` (GO Biological Process): describes the larger cellular or
  physiological role played by a gene product, often in coordination
  with other genes

- `mf` (GO Molecular Function): describes activities performed by gene
  products at the molecular level, such as binding or catalysis

- `h` (MSigDB Hallmark): summarizes well-defined biological states

- `c1` (MSigDB Positional): gene sets corresponding to human chromosome
  cytogenetic bands

- `c2` (MSigDB Curated): selected gene sets from online pathway
  databases and biomedical literature (KEGG and BioCarta removed for
  licensing reasons)

- `c3` (MSigDB Regulatory Target): potential targets of regulation by
  transcription factors or miRNAs

- `c4` (MSigDB Computational): defined by mining cancer-oriented
  expression data (KEGG pathways removed for licensing reasons)

- `c6` (MSigDB Oncogenic Signature): signatures of cellular pathways
  often dysregulated in cancer

- `c7` (MSigDB Immunologic Signature): gene sets associated with immune
  disease

- `c8` (MSigDB Cell Type Signature): curated markers for cell types
  identified in single-cell sequencing studies of human tissue

## Details

[Gene Ontology](https://geneontology.org/) data (retrieved via the
`GO.db` R package, v3.23.1) is made available under the terms of the [CC
BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).

[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) data (version
2026.1.Hs, retrieved via the `msigdbr` R package v26.1.0) is copyright ©
2004–2025 Broad Institute, Inc., Massachusetts Institute of Technology,
and Regents of the University of California, and is made available under
the terms of the [CC BY 4.0
license](https://creativecommons.org/licenses/by/4.0/).

Data was retrieved from these resources on 06/01/2026.

## Author

Amanda Hiser

## Examples

``` r
head(pathway_map)
#>   gene_symbol entrez_id pathway_id group_code
#> 1        PIGV     55650 GO:0000009         mf
#> 2       ALG12     79087 GO:0000009         mf
#> 3        ALG2     85365 GO:0000009         mf
#> 4       PARP1       142 GO:0000012         bp
#> 5       ERCC8      1161 GO:0000012         bp
#> 6       ERCC6      2074 GO:0000012         bp
length(unique(pathway_map$pathway_id))
#> [1] 36791
table(pathway_map$group_code)
#> 
#>      bp      c1      c2      c3      c4      c6      c7      c8 
#> 1144836   44167  606696 1178661   49983   30747  993070  157637 
#>       h      mf 
#>    7331  271232 
```
