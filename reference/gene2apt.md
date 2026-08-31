# Convert Gene Identifiers to SomaScan Analytes (v2)

Takes a set of gene identifiers and returns associated SomaScan
analytes, given a `col_meta` object returned from
[`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html).

## Usage

``` r
gene2apt(
  x,
  col_meta_df = NULL,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  collapse = TRUE,
  sep = "|"
)
```

## Arguments

- x:

  Character. A vector of gene identifiers.

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

  Character. Type of gene identifier provided in `x`. Options are
  "EntrezGeneSymbol" or "EntrezGeneID". Default is "EntrezGeneSymbol".

- collapse:

  Logical. Should values be collapsed to maintain the length of the
  original vector? Default is TRUE. Most useful when the results will be
  used as a data frame column. When FALSE, all genes associated with a
  given `AptName` will be returned as individual elements in the output
  vector. See examples.

- sep:

  Character used to separate collapsed identifiers. Default is "\|".

## Value

A character vector of SomaScan identifiers (in `AptName` format).

## Author

Amanda Hiser, Alex Poole

## Examples

``` r
anno <- SomaDataIO::getAnalyteInfo(example_data_11k)
genes <- withr::with_seed(123, sample(anno$EntrezGeneSymbol, 5))
gene2apt(x = genes, col_meta_df = anno)
#> [1] "seq.18483.36" "seq.18898.36" "seq.17766.5"  "seq.11383.41"
#> [5] "seq.24648.8" 

# By default, one string for each gene will be returned.
# Genes targeted by >1 analytes will have analytes collapsed via `sep`
gene2apt(x = "IL12B", col_meta_df = anno)
#> [1] "seq.10365.132|seq.10367.62|seq.13733.5"

# Use `collapse = FALSE` to return all values individually
gene2apt(x = "IL12B", col_meta_df = anno, collapse = FALSE)
#> [1] "seq.10365.132" "seq.10367.62"  "seq.13733.5"  
```
