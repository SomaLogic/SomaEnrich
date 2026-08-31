# Convert SomaScan Analytes to Genes

Transform a SomaScan analyte into the gene associated with the analyte's
target.

## Usage

``` r
apt2gene(
  x,
  col_meta_df = NULL,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  collapse = TRUE,
  verbose = interactive()
)
```

## Arguments

- x:

  Character. A vector of SomaScan analyte identifiers in `AptName`
  (seq.1234.56) format. Blank ("") values will be coerced to NA. For
  more information about SomaScan identifiers and their formats, please
  see
  [SomaDataIO::SeqId](https://somalogic.github.io/SomaDataIO/reference/SeqId.html).

- col_meta_df:

  Data frame of column metadata, likely created by
  [`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html).
  This object should contain SomaScan menu information that can be used
  to map between `AptNames` and Entrez gene identifiers. At minimum, it
  must contain columns named "AptName" (containing SomaScan identifiers
  in `AptName` format) and either "EntrezGeneSymbol" or "EntrezGeneID".

  If not specified, column metadata will be retrieved from
  `example_data_11k` by default.

- id_type:

  Character. Type of gene identifier to be used in the output. Options
  are "EntrezGeneSymbol" or "EntrezGeneID". Default is
  "EntrezGeneSymbol".

- collapse:

  Logical. Should values be collapsed to maintain the length of the
  original vector? Default is TRUE. Most useful when the results will be
  used as a data frame column. When FALSE, all genes associated with a
  given `AptName` will be returned as individual elements in the output
  vector. See examples.

- verbose:

  Logical. Should function messages be printed to the console?

## Value

A character vector of genes associated with the values provided to `x`.
`NA` is returned when a match is not found.

## See also

[`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html)

## Author

Amanda Hiser, Alex Poole

## Examples

``` r
anno <- SomaDataIO::getAnalyteInfo(SomaDataIO::example_data)
apts <- withr::with_seed(123,
                         sample(SomaDataIO::getAnalytes(SomaDataIO::example_data), 10))
apts <- c(apts, "seq.10367.62") # Adding IL12
apt2gene(x = apts, col_meta_df = anno)
#>  [1] "SRGN"        "CHEK2"       "SH3GL3"      "PHF3"       
#>  [5] "CCDC90B"     "ECM1"        "PDE5A"       "FAM159A"    
#>  [9] "SMPDL3A"     "MMRN2"       "IL12A IL12B"

# Adding results to a data frame
df <- withr::with_seed(321, data.frame(AptName = apts,
                                       Value = rnorm(11)))
df$Gene <- apt2gene(x = df$AptName, col_meta_df = anno)
df
#>         AptName      Value        Gene
#> 1  seq.19251.56  1.7049032        SRGN
#> 2  seq.19328.51 -0.7120386       CHEK2
#> 3  seq.18318.98 -0.2779849      SH3GL3
#> 4  seq.11544.39 -0.1196490        PHF3
#> 5   seq.7792.58 -0.1239606     CCDC90B
#> 6   seq.3366.51  0.2681838        ECM1
#> 7   seq.16805.5  0.7268415       PDE5A
#> 8  seq.13431.74  0.2331354     FAM159A
#> 9   seq.4771.10  0.3391139     SMPDL3A
#> 10 seq.9723.105 -0.5519147       MMRN2
#> 11 seq.10367.62  0.3477014 IL12A IL12B

# Using `collapse = FALSE` to separate all returned values
apt2gene(x = df$AptName, col_meta_df = anno, collapse = FALSE)
#>  [1] "SRGN"    "CHEK2"   "SH3GL3"  "PHF3"    "CCDC90B" "ECM1"   
#>  [7] "PDE5A"   "FAM159A" "SMPDL3A" "MMRN2"   "IL12A"   "IL12B"  
```
