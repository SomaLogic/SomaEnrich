# Example Aptamer-Based Pathway

Example biological pathway comprised of SomaScan analyte identifiers
(`AptNames`). Pathway is a subset of "HALLMARK_INFLAMMATORY_RESPONSE"
(M5932) from MSigDB, converted to `AptNames` via
[`gene2apt()`](https://somalogic.github.io/SomaEnrich/reference/gene2apt.md).
Only genes mapping to an analyte in the SomaScan menu were retained.

## Source

[gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INFLAMMATORY_RESPONSE](https://somalogic.github.io/SomaEnrich/reference/gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INFLAMMATORY_RESPONSE)

## Examples

``` r
head(ex_apt_pathway, n = 20)
#>  [1] "seq.9843.5"   "seq.3174.2"   "seq.9057.19"  "seq.15627.83"
#>  [5] "seq.2867.52"  "seq.9835.16"  "seq.9312.8"   "seq.10574.10"
#>  [9] "seq.3485.28"  "seq.13381.49" "seq.10550.37" "seq.16919.1" 
#> [13] "seq.11214.40" "seq.22128.8"  "seq.11494.4"  "seq.21577.35"
#> [17] "seq.12650.43" "seq.12420.10" "seq.19273.3"  "seq.5230.99" 

length(ex_apt_pathway)
#> [1] 65
```
