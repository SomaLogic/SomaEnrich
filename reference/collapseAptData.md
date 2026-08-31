# Transform Aptamer-Centric Data to Gene-Centric Data

Utility function to convert a `SeqId` or `AptName`-centric data frame,
with one row per analyte, to a gene-centric data frame, with one row per
gene. In cases where more than one analyte is associated with a given
gene, the analytes will be collapsed into one row, with a delimiter to
separate unique values.

## Usage

``` r
collapseAptData(
  df,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  sep = "|",
  apt_order = NULL
)
```

## Arguments

- df:

  Data frame to be converted, most often the output of
  [`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html).
  Must contain, at minimum, an "AptName" column, and the column
  specified by `id_type`.

- id_type:

  Character. Type of gene identifier to be used in the output. Options
  are "EntrezGeneSymbol" or "EntrezGeneID". Default is
  "EntrezGeneSymbol".

- sep:

  Character used to separate collapsed SomaScan identifiers. Default is
  "\|". If NULL, values will *not* be collapsed, and a character vector
  will be returned in each column field when \>1 SomaScan identifier is
  associated with a given gene. See examples.

- apt_order:

  Character. The name of a column in `df` containing **numeric**
  ordering information for collapsed SomaScan identifiers. By default,
  no ordering is used. Should only be specified when multiple aptamers
  are associated with a given gene, and a particular order is desired
  when the values are pasted together into a single, delimited string.

## Value

A data frame with one row per unique gene.

## Author

Amanda Hiser, Alex Poole

## Examples

``` r
df <- SomaDataIO::getAnalyteInfo(example_data_11k)
res <- collapseAptData(df = df)

# Pulling a couple example genes
res[res$EntrezGeneSymbol %in% c("NOTCH1", "APOE", "CAMP"), ]
#>      EntrezGeneSymbol                                         AptName
#> 310              APOE seq.2418.55|seq.2937.10|seq.2938.55|seq.5312.49
#> 750              CAMP                        seq.15481.45|seq.9384.17
#> 3889           NOTCH1                                      seq.5107.7

# Specifying a different delimiter
res_delim <- collapseAptData(df, sep = ", ")
res_delim[res_delim$EntrezGeneSymbol %in% c("NOTCH1", "APOE", "CAMP"), ]
#>      EntrezGeneSymbol
#> 310              APOE
#> 750              CAMP
#> 3889           NOTCH1
#>                                                 AptName
#> 310  seq.2418.55, seq.2937.10, seq.2938.55, seq.5312.49
#> 750                           seq.15481.45, seq.9384.17
#> 3889                                         seq.5107.7

# Using a ranking metric to reorder the analytes
df$Analyte_Rank <- sample(1:nrow(df), nrow(df), replace = FALSE)
res_ord <- collapseAptData(df, apt_order = "Analyte_Rank")
res_ord[res_delim$EntrezGeneSymbol %in% c("NOTCH1", "APOE", "CAMP"), ]
#>      EntrezGeneSymbol                                         AptName
#> 310              APOE seq.2418.55|seq.2938.55|seq.2937.10|seq.5312.49
#> 750              CAMP                        seq.15481.45|seq.9384.17
#> 3889           NOTCH1                                      seq.5107.7

#' Setting `sep = NULL` to create a list-column of SomaScan IDs
res_list <- collapseAptData(df, apt_order = "Analyte_Rank", sep = NULL)
res_list$AptName["APOE"]
#> $APOE
#> [1] "seq.2418.55" "seq.2938.55" "seq.2937.10" "seq.5312.49"
#> 
```
