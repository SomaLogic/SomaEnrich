# Example Gene-Based Pathway

Example biological pathway comprised of gene identifiers (Entrez Gene
IDs). Pathway is a subset of "HALLMARK_INFLAMMATORY_RESPONSE" (M5932)
from MSigDB.

## Source

[gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INFLAMMATORY_RESPONSE](https://somalogic.github.io/SomaEnrich/reference/gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INFLAMMATORY_RESPONSE)

## Examples

``` r
head(ex_gene_pathway, n = 20)
#>  [1] "ABCC4"  "ABHD2"  "ACSL3"  "AKAP12" "APPBP2" "ARID5B" "AZGP1" 
#>  [8] "CCND3"  "CDC14B" "CDK6"   "CENPN"  "DBI"    "DHCR24" "ELK4"  
#> [15] "FKBP5"  "GNAI3"  "HERC3"  "HMGCR"  "HMGCS1" "KLK2"  

length(ex_gene_pathway)
#> [1] 43
```
