# Convert SomaScan ID-Based Feature Set into Gene Set

Utility function to convert a set of SomaScan identifiers (i.e.
`AptNames`) into genes, while also filtering to avoid potential gene
duplication produced by multiple SomaScan analytes mapping to the same
gene.

Mapping between genes and SomaScan analytes is not always 1:1, and it's
possible for multiple analytes to map to the same gene. For example,
`13699-6`, `4330-4`, and `8468-19` all map to prostate-specific antigen
(PSA). If all of these analytes were included in a feature set, PSA
would be represented 3 times in that set. However, when performing
over-representation analysis, pathway analysis, or GSEA, it's crucial to
avoid instances where a single feature is represented multiple times in
the same feature set. This can lead to artificial over-representation
when the duplicated feature is counted multiple times.

After converting SomaScan identifiers to their associated gene, this
function will ensure that each gene is represented only *once* in the
gene set or pathway. However, because 1:many mapping may be present
between identifiers, the pathway returned by `aptPath2genePath()` may be
longer than the input aptamer list. Because of this feature, it is
recommended to use this function instead of
[`apt2gene()`](https://somalogic.github.io/SomaEnrich/reference/apt2gene.md)
when converting entire gene sets or pathways. See the examples for more
details.

## Usage

``` r
aptPath2genePath(
  path_list,
  col_meta_df,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  verbose = interactive()
)
```

## Arguments

- path_list:

  Named list containing character vectors of SomaScan `AptNames`.

- col_meta_df:

  Data frame of column metadata, likely created by
  [`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html).
  This object should contain SomaScan menu information that can be used
  to map between `AptNames` and Entrez gene identifiers. At minimum, it
  must contain columns named "AptName" (containing SomaScan identifiers
  in `AptName` format) and either "EntrezGeneSymbol" or "EntrezGeneID".

- id_type:

  Character. Desired gene identifier to be used in the output. Options
  include "EntrezGeneSymbol" or "EntrezGeneID". Default is
  "EntrezGeneSymbol".

- verbose:

  Logical. Should progress messages be printed to the console?

## Value

The input `pathway` vector, modified to use gene identifiers instead of
`AptNames`. `AptNames` that do not have an associated gene in the
SomaScan menu (determined using the column metadata in the `soma_adat`
object provided to `adat =`) will be dropped.

## Author

Amanda Hiser

## Examples

``` r
# Gene symbols are returned by default
aptPath2genePath(path_list   = list(Example = ex_apt_pathway),
                 col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k))
#> $Example
#>  [1] "SPDEF"    "BMPR1B"   "B2M"      "INSIG1"   "DNAJB9"   "PTK2B"   
#>  [7] "ELL2"     "GPD1L"    "TNFAIP8"  "GNAI3"    "HOMER2"   "B4GALT1" 
#> [13] "HMGCS1"   "UAP1"     "KLK3"     "HSD17B14" "SORD"     "KRT19"   
#> [19] "AKT1"     "LMAN1"    "RAB4A"    "DBI"      "IDI1"     "SMS"     
#> [25] "SRP19"    "NDRG1"    "GSR"      "MYL12A"   "MERTK"    "ITGAV"   
#> [31] "ITGB3"    "ITGB6"    "ITGB8"    "ITGB1"    "SGK1"     "KLK2"    
#> [37] "FKBP5"    "ELK4"     "PDLIM5"   "INPP4B"   "ZBTB10"   "XRCC6"   
#> [43] "UBE2I"    "ADAMTS1"  "RPS6KA3"  "PA2G4"    "ITGB5"    "HPGD"    
#> [49] "HMGCR"    "LIFR"     "UBE2J1"   "PMEPA1"   "VAPA"     "ADRM1"   
#> [55] "AZGP1"    "ALDH1A3"  "ACTN1"    "SAT1"    
#> 

# When a SomaScan assay analyte is associated with more than one gene,
# all genes will be returned
aptPath2genePath(path_list = list(Single_Analyte = "seq.16927.9"),
                 col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k))
#> $Single_Analyte
#> [1] "F13A1" "F13B" 
#> 
```
