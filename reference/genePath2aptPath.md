# Convert Gene Set into SomaScan ID-Based Feature Set

Converts gene pathways to SomaScan identifiers. All analytes that map to
a given gene are retained; if multiple analytes map to the same gene,
all will be present in the output.

## Usage

``` r
genePath2aptPath(
  path_list,
  col_meta_df,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  verbose = interactive()
)
```

## Arguments

- path_list:

  Named list of character vectors containing genes.

- col_meta_df:

  Data frame of column metadata, likely created by
  [`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html).
  This object should contain SomaScan menu information that can be used
  to map between `AptNames` and Entrez gene identifiers. At minimum, it
  must contain columns named "AptName" (containing SomaScan identifiers
  in `AptName` format) and either "EntrezGeneSymbol" or "EntrezGeneID".

- id_type:

  Character. Type of gene identifier used in `path_list`. Options are
  "EntrezGeneSymbol" or "EntrezGeneID". Default is "EntrezGeneSymbol".

- verbose:

  Logical. Show status messages and (if applicable) a progress bar?

## Value

List of named character vectors (names = genes, values = AptNames). When
multiple analytes map to the same gene, all are retained as separate
elements with the same name.

## Author

Amanda Hiser

## Examples

``` r
genePath2aptPath(path_list = list(Example = ex_gene_pathway), 
                col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k))
#> $Example
#>         TNFAIP8           GNAI3          HMGCS1            UAP1 
#>   "seq.12563.2"  "seq.12650.43"  "seq.13496.19"   "seq.13580.2" 
#>            SORD             DBI            SGK1            KLK2 
#>  "seq.15447.45"   "seq.16919.1" "seq.21144.160"  "seq.21444.40" 
#>           FKBP5            ELK4          PDLIM5          ZBTB10 
#>  "seq.21577.35"   "seq.22128.8"   "seq.23703.8"   "seq.25075.2" 
#>           XRCC6         RPS6KA3           HMGCR          UBE2J1 
#>    "seq.2835.1"   "seq.3469.74"   "seq.5230.99"   "seq.6900.30" 
#>           AZGP1 
#>    "seq.9312.8" 
#> 
```
