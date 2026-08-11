# `somaORA()` produces the expected output with 'H' resource

    Code
      somaORA(features = deg, universe = bg, resource = "h", verbose = TRUE)
    Output
      -- Pre-Processing --------------------------------------------------------------
    Message
      x Duplicate elements found in `universe`. 
      > These elements will be removed to create a unique vector.
    Output
      -- Performing Enrichment Analysis ----------------------------------------------
    Message
      > Running ORA...
      v Done!
    Output
         resource_code pathway_id                                    pathway
      1              h      M5916                    HALLMARK_APICAL_SURFACE
      2              h      M5892           HALLMARK_CHOLESTEROL_HOMEOSTASIS
      3              h      M5905                      HALLMARK_ADIPOGENESIS
      4              h      M5902                         HALLMARK_APOPTOSIS
      5              h      M5895        HALLMARK_WNT_BETA_CATENIN_SIGNALING
      6              h      M5919                HALLMARK_HEDGEHOG_SIGNALING
      7              h      M5938   HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
      8              h      M5950               HALLMARK_ALLOGRAFT_REJECTION
      9              h      M5949                        HALLMARK_PEROXISOME
      10             h      M5947               HALLMARK_IL2_STAT5_SIGNALING
      11             h      M5932             HALLMARK_INFLAMMATORY_RESPONSE
      12             h      M5953                 HALLMARK_KRAS_SIGNALING_UP
      13             h      M5915                   HALLMARK_APICAL_JUNCTION
      14             h      M5897           HALLMARK_IL6_JAK_STAT3_SIGNALING
      15             h      M5930 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
      16             h      M5956                 HALLMARK_KRAS_SIGNALING_DN
      17             h      M5941                    HALLMARK_UV_RESPONSE_UP
      18             h      M5909                        HALLMARK_MYOGENESIS
      19             h      M5913         HALLMARK_INTERFERON_GAMMA_RESPONSE
      20             h      M5891                           HALLMARK_HYPOXIA
      21             h      M5934             HALLMARK_XENOBIOTIC_METABOLISM
      22             h      M5890           HALLMARK_TNFA_SIGNALING_VIA_NFKB
      23             h      M5937                        HALLMARK_GLYCOLYSIS
      24             h      M5921                        HALLMARK_COMPLEMENT
      25             h      M5908                 HALLMARK_ANDROGEN_RESPONSE
      26             h      M5944                      HALLMARK_ANGIOGENESIS
      27             h      M5948              HALLMARK_BILE_ACID_METABOLISM
      28             h      M5946                       HALLMARK_COAGULATION
      29             h      M5898                        HALLMARK_DNA_REPAIR
      30             h      M5925                       HALLMARK_E2F_TARGETS
      31             h      M5906           HALLMARK_ESTROGEN_RESPONSE_EARLY
      32             h      M5907            HALLMARK_ESTROGEN_RESPONSE_LATE
      33             h      M5935             HALLMARK_FATTY_ACID_METABOLISM
      34             h      M5901                    HALLMARK_G2M_CHECKPOINT
      35             h      M5945                   HALLMARK_HEME_METABOLISM
      36             h      M5911         HALLMARK_INTERFERON_ALPHA_RESPONSE
      37             h      M5893                   HALLMARK_MITOTIC_SPINDLE
      38             h      M5924                  HALLMARK_MTORC1_SIGNALING
      39             h      M5926                    HALLMARK_MYC_TARGETS_V1
      40             h      M5928                    HALLMARK_MYC_TARGETS_V2
      41             h      M5903                   HALLMARK_NOTCH_SIGNALING
      42             h      M5936         HALLMARK_OXIDATIVE_PHOSPHORYLATION
      43             h      M5939                       HALLMARK_P53_PATHWAY
      44             h      M5957               HALLMARK_PANCREAS_BETA_CELLS
      45             h      M5923           HALLMARK_PI3K_AKT_MTOR_SIGNALING
      46             h      M5910                 HALLMARK_PROTEIN_SECRETION
      47             h      M5951                   HALLMARK_SPERMATOGENESIS
      48             h      M5896                HALLMARK_TGF_BETA_SIGNALING
      49             h      M5922         HALLMARK_UNFOLDED_PROTEIN_RESPONSE
      50             h      M5942                    HALLMARK_UV_RESPONSE_DN
               pval      padj foldEnrichment overlap size    overlapFeatures
      1  0.02361826 0.8838912      8.3337679       2   26         GHRL, THY1
      2  0.05458940 0.8838912      5.2848284       2   41      ALCAM, LGALS3
      3  0.06752457 0.8838912      3.1864407       3  102 ADIPOQ, ENPP2, LEP
      4  0.07071130 0.8838912      3.1251630       3  104  CD2, LGALS3, SOD2
      5  0.18482035 1.0000000      4.9244992       1   22              WNT5B
      6  0.20726568 1.0000000      4.3335593       1   25               THY1
      7  0.25034905 1.0000000      3.4948059       1   31               SOD2
      8  0.14492525 1.0000000      2.2728458       3  143     CD2, LTB, THY1
      9  0.37214893 1.0000000      2.1667797       1   50               SOD2
      10 0.27335625 1.0000000      1.9520537       2  111         ALCAM, LTB
      11 0.29742022 1.0000000      1.8362540       2  118          GPC3, LTA
      12 0.30428054 1.0000000      1.8056497       2  120         IGF2, RELN
      13 0.31113042 1.0000000      1.7760489       2  122        CNTN1, THY1
      14 0.50783757 1.0000000      1.4255129       1   76                LTB
      15 0.41795490 1.0000000      1.4069998       2  154          NTM, THY1
      16 0.57253479 1.0000000      1.1905383       1   91               FSHB
      17 0.57653740 1.0000000      1.1775976       1   92               SOD2
      18 0.64259505 1.0000000      0.9848998       1  110              FABP3
      19 0.66545089 1.0000000      0.9259742       1  117               SOD2
      20 0.68089054 1.0000000      0.8880245       1  122               GPC3
      21 0.68389332 1.0000000      0.8808047       1  123              FETUB
      22 0.68981580 1.0000000      0.8667119       1  125               SOD2
      23 0.70969807 1.0000000      0.8207499       1  132               GPC3
      24 0.75291535 1.0000000      0.7271073       1  149             LGALS3
      25 1.00000000 1.0000000      0.0000000       0   53                   
      26 1.00000000 1.0000000      0.0000000       0   29                   
      27 1.00000000 1.0000000      0.0000000       0   41                   
      28 1.00000000 1.0000000      0.0000000       0  114                   
      29 1.00000000 1.0000000      0.0000000       0   70                   
      30 1.00000000 1.0000000      0.0000000       0   96                   
      31 1.00000000 1.0000000      0.0000000       0   80                   
      32 1.00000000 1.0000000      0.0000000       0   96                   
      33 1.00000000 1.0000000      0.0000000       0  101                   
      34 1.00000000 1.0000000      0.0000000       0   80                   
      35 1.00000000 1.0000000      0.0000000       0   84                   
      36 1.00000000 1.0000000      0.0000000       0   50                   
      37 1.00000000 1.0000000      0.0000000       0   66                   
      38 1.00000000 1.0000000      0.0000000       0  111                   
      39 1.00000000 1.0000000      0.0000000       0  119                   
      40 1.00000000 1.0000000      0.0000000       0   21                   
      41 1.00000000 1.0000000      0.0000000       0   19                   
      42 1.00000000 1.0000000      0.0000000       0   78                   
      43 1.00000000 1.0000000      0.0000000       0  101                   
      44 1.00000000 1.0000000      0.0000000       0   22                   
      45 1.00000000 1.0000000      0.0000000       0   78                   
      46 1.00000000 1.0000000      0.0000000       0   52                   
      47 1.00000000 1.0000000      0.0000000       0   56                   
      48 1.00000000 1.0000000      0.0000000       0   32                   
      49 1.00000000 1.0000000      0.0000000       0   57                   
      50 1.00000000 1.0000000      0.0000000       0   75                   

# `somaORA()` produces a progress message when `verbose = TRUE`

    Code
      t <- somaORA(features = deg, universe = bg, resource = "h", verbose = TRUE)
    Output
      -- Pre-Processing --------------------------------------------------------------
    Message
      x Duplicate elements found in `universe`. 
      > These elements will be removed to create a unique vector.
    Output
      -- Performing Enrichment Analysis ----------------------------------------------
    Message
      > Running ORA...
      v Done!

