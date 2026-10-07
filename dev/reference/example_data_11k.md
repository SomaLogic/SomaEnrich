# Example 11K Data Set

The `example_data_11k` object is intended to provide existing and
prospective SomaLogic customers with example data to enable analysis
preparation prior to receipt of SomaScan data, and also for those
generally curious about the SomaScan data deliverable. It is **not**
intended to be used as a control group for studies or provide any
metrics for SomaScan data in general.

## Source

<https://github.com/SomaLogic/SomaLogic-Data/blob/main/example_data_v5.0_plasma.adat>

## Examples

``` r
# S3 print method
example_data_11k
#> ══ SomaScan Data ══════════════════════════════════════════════════════
#>      SomaScan version     v4.1 (7k)
#>      Signal Space         7k
#>      Attributes intact    ✓
#>      Rows                 207
#>      Columns              7632
#>      Clinical Data        36
#>      Features             7596
#> ── Column Meta ────────────────────────────────────────────────────────
#> ℹ SeqId, SeqIdVersion, SomaId, TargetFullName, Target,
#> ℹ UniProt, EntrezGeneID, EntrezGeneSymbol, Organism, Units,
#> ℹ Type, Dilution, PlateScale_Reference, CalReference,
#> ℹ Cal_SS.000005_Set001, ColCheck,
#> ℹ CalQcRatio_SS.000005_Set001_200170, QcReference_200170,
#> ℹ Cal_SS.000005_Set002, CalQcRatio_SS.000005_Set002_200170,
#> ℹ Cal_SS.000005_Set003, CalQcRatio_SS.000005_Set003_200170,
#> ℹ Cal_SS.000005_Set004, CalQcRatio_SS.000005_Set004_200170,
#> ℹ Dilution2
#> ── Tibble ─────────────────────────────────────────────────────────────
#> # A tibble: 207 × 7,633
#>    row_names      PlateId  PlateRunDate ScannerID PlatePosition SlideId
#>    <chr>          <chr>    <chr>        <chr>     <chr>           <dbl>
#>  1 258495800070_8 SS-0000… 2021-11-30   SG152144… A1            2.58e11
#>  2 258495800075_6 SS-0000… 2021-11-30   SG152144… A11           2.58e11
#>  3 258495800068_8 SS-0000… 2021-11-30   SG152144… A12           2.58e11
#>  4 258495800068_7 SS-0000… 2021-11-30   SG152144… A2            2.58e11
#>  5 258495800064_4 SS-0000… 2021-11-30   SG152144… A4            2.58e11
#>  6 258495800072_8 SS-0000… 2021-11-30   SG152144… A5            2.58e11
#>  7 258495800072_5 SS-0000… 2021-11-30   SG152144… A6            2.58e11
#>  8 258495800074_5 SS-0000… 2021-11-30   SG152144… A8            2.58e11
#>  9 258495800067_5 SS-0000… 2021-11-30   SG152144… A9            2.58e11
#> 10 258495800071_1 SS-0000… 2021-11-30   SG152144… B1            2.58e11
#> # ℹ 197 more rows
#> # ℹ 7,627 more variables: Subarray <dbl>, SampleId <chr>,
#> #   SampleType <chr>, PercentDilution <int>, SampleMatrix <chr>,
#> #   Barcode <lgl>, Barcode2d <chr>, SampleName <lgl>,
#> #   SampleNotes <lgl>, AliquotingNotes <lgl>, …
#> ═══════════════════════════════════════════════════════════════════════

# Print header info
print(example_data_11k, show_header = TRUE)
#> ══ SomaScan Data ══════════════════════════════════════════════════════
#>      SomaScan version     v4.1 (7k)
#>      Signal Space         7k
#>      Attributes intact    ✓
#>      Rows                 207
#>      Columns              7632
#>      Clinical Data        36
#>      Features             7596
#> ── Column Meta ────────────────────────────────────────────────────────
#> ℹ SeqId, SeqIdVersion, SomaId, TargetFullName, Target,
#> ℹ UniProt, EntrezGeneID, EntrezGeneSymbol, Organism, Units,
#> ℹ Type, Dilution, PlateScale_Reference, CalReference,
#> ℹ Cal_SS.000005_Set001, ColCheck,
#> ℹ CalQcRatio_SS.000005_Set001_200170, QcReference_200170,
#> ℹ Cal_SS.000005_Set002, CalQcRatio_SS.000005_Set002_200170,
#> ℹ Cal_SS.000005_Set003, CalQcRatio_SS.000005_Set003_200170,
#> ℹ Cal_SS.000005_Set004, CalQcRatio_SS.000005_Set004_200170,
#> ℹ Dilution2
#> ── Header Data ────────────────────────────────────────────────────────
#> # A tibble: 58 × 2
#>    Key                  Value                                          
#>    <chr>                <chr>                                          
#>  1 AdatId               GID-baeca5fc-df4b-4dde-a96c-fa5a80433c23       
#>  2 Version              1.2                                            
#>  3 AssayType            PharmaServices                                 
#>  4 AssayVersion         v4.1                                           
#>  5 AssayRobot           Fluent                                         
#>  6 Legal                Experiment details and data have been processe…
#>  7 CreatedBy            PharmaServices                                 
#>  8 CreatedDate          2022-01-21                                     
#>  9 EnteredBy            Technician1                                    
#> 10 ExpDate              2021-11-30, 2021-12-01                         
#> 11 GeneratedBy          Px (Build:  : ), Canopy_0.1.1, Canopy_0.3      
#> 12 ProcessSteps         Raw RFU, Hyb Normalization, medNormInt (Sample…
#> 13 ProteinEffectiveDate 2021-07-01                                     
#> 14 RunNotes             2 columns ('Age' and 'Sex') have been added to…
#> 15 StudyMatrix          EDTA Plasma                                    
#> # ℹ 43 more rows
#> ═══════════════════════════════════════════════════════════════════════

# View object class
class(example_data_11k)
#> [1] "soma_adat"  "data.frame"

# Retrieve annotations
SomaDataIO::getAnalyteInfo(example_data_11k)
#> # A tibble: 7,596 × 26
#>    AptName      SeqId SeqIdVersion SomaId TargetFullName Target UniProt
#>    <chr>        <chr>        <dbl> <chr>  <chr>          <chr>  <chr>  
#>  1 seq.10000.28 1000…            3 SL019… Beta-crystall… CRBB2  P43320 
#>  2 seq.10001.7  1000…            3 SL002… RAF proto-onc… c-Raf  P04049 
#>  3 seq.10003.15 1000…            3 SL019… Zinc finger p… ZNF41  P51814 
#>  4 seq.10006.25 1000…            3 SL019… ETS domain-co… ELK1   P19419 
#>  5 seq.10008.43 1000…            3 SL019… Guanylyl cycl… GUC1A  P43080 
#>  6 seq.10010.10 1001…            3 SL014… Beclin-1       BECN1  Q14457 
#>  7 seq.10011.65 1001…            3 SL019… Inositol poly… OCRL   Q01968 
#>  8 seq.10012.5  1001…            3 SL014… SAM pointed d… SPDEF  O95238 
#>  9 seq.10013.34 1001…            3 SL025… Fc_MOUSE       Fc_MO… Q99LC4 
#> 10 seq.10014.31 1001…            3 SL007… Zinc finger p… SLUG   O43623 
#> # ℹ 7,586 more rows
#> # ℹ 19 more variables: EntrezGeneID <chr>, EntrezGeneSymbol <chr>,
#> #   Organism <chr>, Units <chr>, Type <chr>, Dilution <chr>,
#> #   PlateScale_Reference <dbl>, CalReference <dbl>,
#> #   Cal_SS.000005_Set001 <dbl>, ColCheck <chr>,
#> #   CalQcRatio_SS.000005_Set001_200170 <dbl>,
#> #   QcReference_200170 <dbl>, Cal_SS.000005_Set002 <dbl>, …
```
