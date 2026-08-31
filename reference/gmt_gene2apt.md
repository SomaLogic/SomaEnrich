# Convert Gene-Based GMT File to AptNames

Reads in a standard GMT file (using gene identifiers) and converts genes
to SomaScan `AptNames`.

## Usage

``` r
gmt_gene2apt(
  gmt_gene_file,
  col_meta_df,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  write = FALSE,
  gmt_apt_file = NULL,
  verbose = interactive()
)
```

## Arguments

- gmt_gene_file:

  Path to the `.gmt` file to be converted.

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

- write:

  Logical. Should the converted data set be written out to a new GMT
  file? Default is FALSE.

- gmt_apt_file:

  The desired file name for the converted file. Required if
  `write = TRUE`. The file will be written to the current directory,
  unless a path to the final file is provided.

- verbose:

  Logical. Show status messages and (if applicable) a progress bar?

## Value

A named list of gene sets/pathways, converted to SomaScan `AptNames`.

## Examples

``` r
if (FALSE) { # \dontrun{
example_gmt <- system.file("extdata", 
                           "msigdb_h_androgen_response.gmt", 
                            package = "SomaEnrich")
gmt_gene2apt(example_gmt, 
             col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k), 
             gmt_apt_file = "msigdb_h_androgen_response_apt.gmt")
} # }
```
