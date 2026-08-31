# Retrieve SomaScan Analytes in a Given GO Term

Retrieve SomaScan Analytes in a Given GO Term

## Usage

``` r
go2apt(x, col_meta_df = NULL, verbose = interactive())
```

## Arguments

- x:

  Character. A single GO term ID. Must contain "GO:" prefix.

- col_meta_df:

  Data frame of column metadata, likely created by
  [`SomaDataIO::getAnalyteInfo()`](https://somalogic.github.io/SomaDataIO/reference/getAnalyteInfo.html).
  This object should contain SomaScan menu information that can be used
  to map between `AptNames` and Entrez gene identifiers. At minimum, it
  must contain columns named "AptName" (containing SomaScan identifiers
  in `AptName` format) and either "EntrezGeneSymbol" or "EntrezGeneID".

  If not specified, column metadata will be retrieved from
  `example_data_11k` by default.

- verbose:

  Logical. Should function messages be printed to the console?

## Value

Character vector of SomaScan analyte identifiers in `AptName` format.
All analytes associated with each gene are returned. For more
information about SomaScan identifiers and their formats, please see
[SomaDataIO::SeqId](https://somalogic.github.io/SomaDataIO/reference/SeqId.html).

## Examples

``` r
go2apt("GO:0019319") # Returns vector of AptNames
#>  [1] "seq.10339.48"  "seq.11083.23"  "seq.11105.171" "seq.11286.78" 
#>  [5] "seq.11327.56"  "seq.11825.27"  "seq.12456.5"   "seq.12649.80" 
#>  [9] "seq.12651.21"  "seq.12960.9"   "seq.13936.24"  "seq.13990.1"  
#> [13] "seq.14006.36"  "seq.14097.86"  "seq.15447.45"  "seq.15524.30" 
#> [17] "seq.15534.26"  "seq.15633.6"   "seq.16606.85"  "seq.16616.137"
#> [21] "seq.17333.20"  "seq.18182.24"  "seq.18185.118" "seq.18235.16" 
#> [25] "seq.19369.17"  "seq.19376.74"  "seq.19797.4"   "seq.19823.75" 
#> [29] "seq.20079.6"   "seq.20128.1"   "seq.21733.11"  "seq.21833.6"  
#> [33] "seq.22055.31"  "seq.22141.59"  "seq.22405.61"  "seq.25041.11" 
#> [37] "seq.3401.8"    "seq.3466.8"    "seq.3554.24"   "seq.3896.5"   
#> [41] "seq.4272.46"   "seq.4309.59"   "seq.4407.10"   "seq.4469.78"  
#> [45] "seq.4883.56"   "seq.4891.50"   "seq.4912.17"   "seq.4963.19"  
#> [49] "seq.5020.50"   "seq.5183.53"   "seq.5245.40"   "seq.7206.20"  
#> [53] "seq.7251.64"   "seq.7831.39"   "seq.9173.21"   "seq.9854.36"  
#> [57] "seq.9867.23"   "seq.9876.20"  

# Apply over a vector of GO terms
sapply(c("GO:2001256", "GO:0000462", "GO:2001259"), go2apt, USE.NAMES = TRUE)
#> $`GO:2001256`
#> [1] "seq.11263.57" "seq.12685.57" "seq.19229.92" "seq.21799.15"
#> [5] "seq.22980.37" "seq.3642.4"   "seq.8073.3"   "seq.8916.32" 
#> [9] "seq.9271.101"
#> 
#> $`GO:0000462`
#> [1] "seq.19166.15" "seq.20388.8"  "seq.22512.18" "seq.24279.30"
#> [5] "seq.25422.58"
#> 
#> $`GO:2001259`
#>  [1] "seq.10785.8"   "seq.11263.57"  "seq.24237.115" "seq.24304.3"  
#>  [5] "seq.25125.7"   "seq.3181.50"   "seq.6495.14"   "seq.6998.106" 
#>  [9] "seq.8916.32"   "seq.9271.101" 
#> 
```
