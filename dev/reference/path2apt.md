# Retrieve SomaScan Analytes in Pathways from `pathway_map`

Takes one or more pathway identifiers from `pathway_map$pathway_id` (any
resource, e.g. GO, MSigDB Hallmark, etc.) and returns the SomaScan
analytes that target the genes in each pathway. All analytes associated
with each gene are returned.

## Usage

``` r
path2apt(
  x,
  col_meta_df = NULL,
  id_type = c("EntrezGeneSymbol", "EntrezGeneID"),
  verbose = interactive()
)
```

## Arguments

- x:

  Character. A vector of pathway identifiers found in
  `pathway_map$pathway_id` (e.g. "GO:0019319", "M5890").

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

  Character. Gene identifier used to match pathway genes to analytes in
  `col_meta_df`; this does not change the output, which is always
  `AptNames`. Options are "EntrezGeneSymbol" or "EntrezGeneID". Entrez
  Gene IDs are robust to gene symbol changes. Default is
  "EntrezGeneSymbol".

- verbose:

  Logical. Should function messages be printed to the console?

## Value

A named list of character vectors, one element per unique pathway in `x`
(in the order provided). Each vector contains unique SomaScan analyte
identifiers in `AptName` format. Pathways with no matching analytes in
`col_meta_df` return `character(0)`. For more information about SomaScan
identifiers and their formats, please see
[SomaDataIO::SeqId](https://somalogic.github.io/SomaDataIO/reference/SeqId.html).

## Author

Amanda Hiser

## Examples

``` r
path2apt("GO:0019319")
#> $`GO:0019319`
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
#> 

# Multiple pathways can be mapped at once
path2apt(c("M41209", "GO:2001256", "M1568"))
#> $M41209
#>  [1] "seq.10085.25" "seq.13556.28" "seq.13741.36" "seq.15668.19"
#>  [5] "seq.22485.1"  "seq.22857.5"  "seq.24701.21" "seq.25087.11"
#>  [9] "seq.2741.22"  "seq.2771.35"  "seq.2871.73"  "seq.3350.53" 
#> [13] "seq.3469.74"  "seq.4249.64"  "seq.6431.68"  "seq.7136.107"
#> [17] "seq.8074.32"  "seq.8086.49" 
#> 
#> $`GO:2001256`
#> [1] "seq.11263.57" "seq.12685.57" "seq.19229.92" "seq.21799.15"
#> [5] "seq.22980.37" "seq.3642.4"   "seq.8073.3"   "seq.8916.32" 
#> [9] "seq.9271.101"
#> 
#> $M1568
#>  [1] "seq.10082.251" "seq.13654.1"   "seq.14158.17"  "seq.15494.11" 
#>  [5] "seq.15540.6"   "seq.15640.54"  "seq.17163.117" "seq.18178.13" 
#>  [9] "seq.18233.10"  "seq.19239.5"   "seq.19259.176" "seq.19264.6"  
#> [13] "seq.19516.10"  "seq.20592.8"   "seq.23903.3"   "seq.25092.32" 
#> [17] "seq.25248.28"  "seq.2625.53"   "seq.3171.57"   "seq.3714.49"  
#> [21] "seq.3800.71"   "seq.4719.58"   "seq.4829.43"   "seq.4970.55"  
#> [25] "seq.5031.10"   "seq.5033.27"   "seq.5532.53"   "seq.7253.6"   
#> [29] "seq.7832.181"  "seq.9507.55"   "seq.9756.6"    "seq.9870.17"  
#> 

# By default, pathway genes are matched to analytes by gene symbol. Matching
# by Entrez Gene ID instead can recover analytes whose gene was renamed
# (e.g. DDX58 -> RIGI)
sym <- path2apt("M5890")
ent <- path2apt("M5890", id_type = "EntrezGeneID")
setdiff(ent[[1L]], sym[[1L]])
#> [1] "seq.12382.2"
```
