# t-test Example Data

Results of performing differential expression (via t-test) with an
example SomaScan dataset (`example_data_11k`). The `Sex` variable was
used to compare males (N = 81) against females (N = 81). The results are
sorted by t-statistic and adjusted p-value.

## Source

<https://somalogic.github.io/SomaDataIO/articles/stat-two-group-comparison.html>

## Details

The SomaDataIO "Two-Group Comparison" workflow was used to produce this
data set, with modifications (different input ADAT, and fold change and
log2 fold change values were added to the final table).

## Examples

``` r
colnames(t_tests)
#>  [1] "AptName"          "SeqId"            "Target"          
#>  [4] "EntrezGeneSymbol" "UniProt"          "formula"         
#>  [7] "t_test"           "t_stat"           "p.value"         
#> [10] "fdr"              "fc"               "log2_fc"         

head(t_tests)
#> # A tibble: 6 × 12
#>   AptName      SeqId  Target EntrezGeneSymbol UniProt formula   t_test 
#>   <chr>        <chr>  <chr>  <chr>            <chr>   <list>    <list> 
#> 1 seq.6580.29  6580-… Pregn… PZP              P20742  <formula> <htest>
#> 2 seq.3032.11  3032-… FSH    CGA|FSHB         P01215… <formula> <htest>
#> 3 seq.16892.23 16892… ENPP2  ENPP2            Q13822  <formula> <htest>
#> 4 seq.2953.31  2953-… Lutei… CGA|LHB          P01215… <formula> <htest>
#> 5 seq.4914.10  4914-… HCG    CGA|CGB3|CGB7    P01215… <formula> <htest>
#> 6 seq.8428.102 8428-… NTRI   NTM              Q9P121  <formula> <htest>
#> # ℹ 5 more variables: t_stat <dbl>, p.value <dbl>, fdr <dbl>,
#> #   fc <dbl>, log2_fc <dbl>

# Formula used to generate results for a given analyte
t_tests$formula[[1]]
#> seq.6580.29 ~ Sex
#> <environment: 0x12f34d5e0>
```
