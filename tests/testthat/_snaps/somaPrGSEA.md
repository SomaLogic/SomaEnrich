# `somaPrGSEA()` has expected results with 'H' resource

    Code
      res$results
    Output
         resource_code pathway_id                                    pathway
      1              h      M5908                 HALLMARK_ANDROGEN_RESPONSE
      2              h      M5945                   HALLMARK_HEME_METABOLISM
      3              h      M5935             HALLMARK_FATTY_ACID_METABOLISM
      4              h      M5926                    HALLMARK_MYC_TARGETS_V1
      5              h      M5915                   HALLMARK_APICAL_JUNCTION
      6              h      M5898                        HALLMARK_DNA_REPAIR
      7              h      M5934             HALLMARK_XENOBIOTIC_METABOLISM
      8              h      M5948              HALLMARK_BILE_ACID_METABOLISM
      9              h      M5949                        HALLMARK_PEROXISOME
      10             h      M5957               HALLMARK_PANCREAS_BETA_CELLS
      11             h      M5930 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
      12             h      M5938   HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
      13             h      M5891                           HALLMARK_HYPOXIA
      14             h      M5925                       HALLMARK_E2F_TARGETS
      15             h      M5924                  HALLMARK_MTORC1_SIGNALING
      16             h      M5892           HALLMARK_CHOLESTEROL_HOMEOSTASIS
      17             h      M5939                       HALLMARK_P53_PATHWAY
      18             h      M5906           HALLMARK_ESTROGEN_RESPONSE_EARLY
      19             h      M5910                 HALLMARK_PROTEIN_SECRETION
      20             h      M5922         HALLMARK_UNFOLDED_PROTEIN_RESPONSE
      21             h      M5905                      HALLMARK_ADIPOGENESIS
      22             h      M5953                 HALLMARK_KRAS_SIGNALING_UP
      23             h      M5897           HALLMARK_IL6_JAK_STAT3_SIGNALING
      24             h      M5921                        HALLMARK_COMPLEMENT
      25             h      M5936         HALLMARK_OXIDATIVE_PHOSPHORYLATION
      26             h      M5919                HALLMARK_HEDGEHOG_SIGNALING
      27             h      M5956                 HALLMARK_KRAS_SIGNALING_DN
      28             h      M5907            HALLMARK_ESTROGEN_RESPONSE_LATE
      29             h      M5916                    HALLMARK_APICAL_SURFACE
      30             h      M5909                        HALLMARK_MYOGENESIS
      31             h      M5946                       HALLMARK_COAGULATION
      32             h      M5901                    HALLMARK_G2M_CHECKPOINT
      33             h      M5932             HALLMARK_INFLAMMATORY_RESPONSE
      34             h      M5923           HALLMARK_PI3K_AKT_MTOR_SIGNALING
      35             h      M5911         HALLMARK_INTERFERON_ALPHA_RESPONSE
      36             h      M5950               HALLMARK_ALLOGRAFT_REJECTION
      37             h      M5944                      HALLMARK_ANGIOGENESIS
      38             h      M5942                    HALLMARK_UV_RESPONSE_DN
      39             h      M5951                   HALLMARK_SPERMATOGENESIS
      40             h      M5947               HALLMARK_IL2_STAT5_SIGNALING
      41             h      M5902                         HALLMARK_APOPTOSIS
      42             h      M5928                    HALLMARK_MYC_TARGETS_V2
      43             h      M5896                HALLMARK_TGF_BETA_SIGNALING
      44             h      M5913         HALLMARK_INTERFERON_GAMMA_RESPONSE
      45             h      M5937                        HALLMARK_GLYCOLYSIS
      46             h      M5941                    HALLMARK_UV_RESPONSE_UP
      47             h      M5890           HALLMARK_TNFA_SIGNALING_VIA_NFKB
      48             h      M5895        HALLMARK_WNT_BETA_CATENIN_SIGNALING
      49             h      M5893                   HALLMARK_MITOTIC_SPINDLE
      50             h      M5903                   HALLMARK_NOTCH_SIGNALING
                pval       padj    log2err         ES        NES
      1  0.005621573 0.09339826 0.40701792 -0.4648078 -1.7368257
      2  0.003702097 0.09339826 0.43170770 -0.4008774 -1.6738785
      3  0.007471861 0.09339826 0.40701792 -0.3528691 -1.5116464
      4  0.004962902 0.09339826 0.40701792 -0.3365589 -1.4915343
      5  0.021898337 0.21898337 0.35248786  0.3404886  1.4015395
      6  0.040672500 0.33893750 0.32177592 -0.3562255 -1.4310200
      7  0.050632797 0.36166284 0.32177592 -0.2969329 -1.3190582
      8  0.125603865 0.44992175 0.19578900 -0.3593385 -1.2851213
      9  0.099750623 0.44992175 0.22496609 -0.3414397 -1.2664179
      10 0.193043478 0.44992175 0.12944289  0.4265008  1.2465416
      11 0.081942337 0.44992175 0.19189224  0.2907501  1.2413036
      12 0.159705160 0.44992175 0.17374784 -0.3712130 -1.2357298
      13 0.090410959 0.44992175 0.24891111 -0.2721964 -1.2143716
      14 0.121621622 0.44992175 0.21140019 -0.2840089 -1.2052739
      15 0.124309392 0.44992175 0.21140019 -0.2738880 -1.1992106
      16 0.193877551 0.44992175 0.12750532  0.3487094  1.1930436
      17 0.137203166 0.44992175 0.19578900 -0.2777818 -1.1921013
      18 0.139896373 0.44992175 0.19189224 -0.2879772 -1.1862058
      19 0.170483461 0.44992175 0.17093234 -0.3164445 -1.1824434
      20 0.162849873 0.44992175 0.17520405 -0.3074389 -1.1789761
      21 0.178515008 0.44992175 0.12814292  0.2940856  1.1754358
      22 0.215962441 0.44992175 0.11426650  0.2757357  1.1362408
      23 0.212276215 0.44992175 0.15214492 -0.2783519 -1.1352570
      24 0.198300283 0.44992175 0.16693385 -0.2451062 -1.1206170
      25 0.242346939 0.48469388 0.14122512 -0.2702661 -1.1080627
      26 0.263157895 0.50607287 0.10672988  0.3910111  1.1721337
      27 0.283413849 0.52484046 0.09889030  0.2805933  1.0974360
      28 0.301587302 0.53854875 0.12750532 -0.2508583 -1.0667270
      29 0.333907057 0.55651176 0.09314546  0.3557079  1.0847135
      30 0.332425068 0.55651176 0.12267919 -0.2422186 -1.0574374
      31 0.365122616 0.58890744 0.11623415 -0.2381198 -1.0395434
      32 0.416230366 0.63400360 0.10512513 -0.2483022 -1.0164458
      33 0.422222222 0.63400360 0.10797236 -0.2294141 -1.0137917
      34 0.431122449 0.63400360 0.10135074 -0.2451615 -1.0015684
      35 0.454242928 0.64891847 0.07455008  0.2830830  1.0058718
      36 0.506172840 0.68401735 0.06538342  0.2313517  0.9724776
      37 0.504878049 0.68401735 0.08943668 -0.2941429 -0.9609664
      38 0.561374795 0.71415492 0.06321912  0.2527322  0.9555005
      39 0.586734694 0.71415492 0.08359906 -0.2490156 -0.9434754
      40 0.614173228 0.71415492 0.05712585  0.2296153  0.9285546
      41 0.611727417 0.71415492 0.05760911  0.2309865  0.9258519
      42 0.593896714 0.71415492 0.07850290 -0.3018949 -0.9189616
      43 0.578151261 0.71415492 0.06307904  0.2840456  0.9162739
      44 0.629120879 0.71491009 0.08383611 -0.2125306 -0.9423343
      45 0.702312139 0.78034682 0.08063885 -0.2033762 -0.9074187
      46 0.754385965 0.81998474 0.04821497  0.2192686  0.8641392
      47 0.815126050 0.84908964 0.07130530 -0.1932648 -0.8585358
      48 0.800000000 0.84908964 0.04949049  0.2598921  0.7595915
      49 0.859504132 0.87704503 0.04388798  0.2131412  0.7912584
      50 0.931818182 0.93181818 0.05479395 -0.2198620 -0.6557584
                                                                                                                                                                                                                                                                                                                                                                                                                    leadingEdge
      1                                                                                                                                                                                                                                                                                                                                                                  KLK3, HMGCS1, MERTK, AZGP1, SORD, ELL2, B4GALT1, PTK2B
      2                                                                                                                                                                                                        ACP5, SDCBP, BPGM, CA1, BLVRB, HAGH, HMBS, PSMD9, MINPP1, ELL2, HDGF, DCUN1D1, NR3C1, TFRC, HEBP1, FTCD, TCEA1, USP15, GCLM, HTATIP2, UBAC1, SLC25A38, CTSB, ERMAP, HTRA2, TOP1, MAP2K3, ALAD, FOXO3, HBZ, EPB41
      3                                                                                                                                                                              HMGCS1, CA6, AUH, GLUL, S100A10, ALDH1A1, APEX1, ADH1C, HSDL2, ENO3, UGDH, CEL, GRHPR, PDHB, CBR3, G0S2, HPGD, GAPDHS, ACOT8, ECI1, GPD1, ALAD, PPARA, YWHAH, UBE2L6, ODC1, ACADVL, MIF, HADH, GCDH, PRDX6, UROS, MDH2, HSPH1, IDH1, ALDOA
      4  NPM1, HNRNPD, CUL1, RACK1, HDGF, SNRPB2, APEX1, GLO1, RPS3, HNRNPA2B1, RPS5, PHB2, SSBP1, ERH, TXNL4A, YWHAE, XRCC6, PCNA, SNRPA, SSB, CCNA2, HDAC2, SRPK1, PA2G4, PSMA7, HNRNPA1, ODC1, TYMS, KARS1, PPM1G, EIF4A1, NME1, EIF4E, NHP2, NAP1L1, CTPS1, COX5A, PSMA2, ACP1, SNRPA1, PSMA1, RANBP1, RAN, PCBP1, EIF4G2, HSPE1, NOP16, EIF2S2, HSPD1, G3BP1, UBA2, NCBP1, KPNB1, PRDX3, UBE2E1, EIF2S1, PPIA, H2AZ1, PGK1
      5                                                                                                                                                                                                CNTN1, THY1, TNFRSF11B, CADM2, VCAN, CDH1, ADAMTS5, CDSN, CD209, CADM3, CD86, VCL, CDH4, SLIT2, ADAM15, NLGN2, NEGR1, CD34, CALB2, GRB7, LAMC2, COL9A1, NLGN3, AMH, SLC30A3, ICAM2, NECTIN1, ICAM5, CDH15, PARVA, CX3CL1
      6                                                                                                                                                                                                                SDCBP, POLB, TK2, CDA, AK1, NUDT9, TMED2, CETN2, TP53, POLR2I, GMPR2, PCNA, RPA2, POLR2C, SRSF6, RALA, ITPA, ERCC4, ERCC1, FEN1, TYMS, MPG, ARL6IP1, GTF2A2, POLR2F, NME1, PNP, GUK1, RRM2B, NME4, EIF1B
      7                                                                                                                                                                                                                   BLVRB, ASL, AHCY, HMOX1, DCXR, TKFC, TPST1, PTGR1, POR, HSD11B1, ALDH2, NQO1, IL1R1, CDA, ADH1C, AKR1C2, UGDH, CES1, APOE, MTHFD1, SAR1B, GSS, GSTT2, FBP1, UPP1, KYNU, HACL1, RBP4, FAS, FMO3, CCL25
      8                                                                                                                                                                                                                                                                                   ALDH1A1, APOA1, SULT2B1, HAO1, PHYH, SOD1, GCLM, HACL1, SLC27A2, PEX26, LCK, CAT, PECR, AGXT, IDH1, AKR1D1, PIPOX, MLYCD, GNMT, NR3C2
      9                                                                                                                                                                                                                                                                                                                        ALDH1A1, SEMA3C, SULT2B1, PRDX1, ALB, SIAH1, HRAS, CEL, SOD1, SLC27A2, ACOT8, ERCC1, YWHAH, FDPS
      10                                                                                                                                                                                                                                                                                                                                                                        IAPP, SST, PCSK1, HNF1A, SCGN, GCG, SYT13, PAX4
      11                                                                                 NTM, THY1, APLP1, IGFBP3, MATN2, LUM, SFRP4, TNFRSF11B, NOTCH2, VCAN, LRRC15, GREM1, PDGFRB, CAPG, THBS2, RGS4, FBLN1, SLIT2, LRP1, PCOLCE2, MMP14, POSTN, INHBA, ID2, PTHLH, IGFBP4, LAMC2, TFPI2, SLIT3, P3H1, JUN, FMOD, LOXL2, COL6A3, NT5E, FAP, GAS1, PLOD3, PLAUR, BGN, CXCL12, LGALS1, MGP, ITGB5, PTX3, IGFBP2, PFN2, TNFAIP3
      12                                                                                                                                                                                                                                                                                                                                                FTL, NQO1, LSP1, PRDX1, GLRX, SOD1, GCLM, OXSR1, CAT, PRDX6, PTPA, GPX3
      13                                                                                    HMOX1, STC1, S100A4, TGM2, TGFBI, ZFP36, INHA, GPC1, NR3C1, SULT2B1, LDHC, ENO3, IER3, ALDOB, GRHPR, GLRX, EGFR, PGF, FBP1, CDKN1B, TGFB3, PKLR, CHST2, SLC25A1, DCN, GAPDHS, ANXA2, FOXO3, GAPDH, TPST2, EFNA3, BTG1, HEXA, MIF, SERPINE1, IL6, GALK1, CCN1, XPNPEP1, BCL2, NDRG1, PDGFB, ADM, JMJD6, ALDOA, PPARGC1A, PGM1, FOSL2
      14                                                                                                                                                                               HNRNPD, PSIP1, PPP1R8, DCK, SMC3, TFRC, HMGA1, POP7, CHEK2, SLBP, TP53, CDKN1B, TIMELESS, XRCC6, PCNA, RPA2, JPT1, EZH2, EXOSC8, HMGB2, ZW10, PA2G4, TMPO, MTHFD2, UNG, PMS2, EED, NME1, BARD1, NAP1L1, CTPS1, AURKB, RANBP1, CCNE1, RAN
      15                                                                                                                                                                                                                                                                                                           HMGCS1, PHGDH, PSPH, STC1, CTSC, HSPA9, CALR, HMBS, PSAT1, SORD, TFRC, CFP, TCEA1, PRDX1, UBE2D3, GLRX, TBK1
      16                                                                                                                                                                                                                                                                                                                                   LGALS3, ALCAM, CLU, ANTXR2, S100A11, ALDOC, ETHE1, SEMA3B, HSD17B7, CD9, LPL, CXCL16
      17                                                                                                                                                                                                            HMOX1, DCXR, APAF1, ANKRA2, S100A10, SERPINB5, S100A4, H1-2, RACK1, SFN, KLK8, AK1, CD82, HRAS, CTSD, IER3, TP53, UPP1, RPS12, PCNA, FAS, FUCA1, TXNIP, ABHD4, STOM, FOXO3, INHBB, TGFB1, BTG1, TPRKB, SAT1
      18                                                                                                                                                                                                              MSMB, BLVRB, CLIC3, ALDH3B1, TGM2, SFN, B4GALT1, PDLIM3, SULT2B1, SH3BP5, IL17RB, CD44, SLC27A2, RETREG1, ASB13, KRT19, IL6ST, GLA, INPP5F, PGR, RHOD, INHBB, MICB, GJA1, MYB, FKBP4, KRT18, BAG1, OLFML3
      19                                                                                                                                                                                                                                                                                                                    CTSC, SNAP23, GNAS, ICA1, TMED2, ARF1, PPT1, EGFR, SOD1, GLA, ZW10, KRT18, RER1, IGF2R, DNM1L, TMX1
      20                                                                                                                                                                                                                                                                                                                                                                     ATF6, NPM1, KHSRP, EDC4, HSPA9, CALR, PSAT1, RPS14
      21                                                                                                                                                                                                                                                                                                        ENPP2, LEP, ADIPOQ, FABP4, OMD, SPARCL1, BCL6, DLD, C3, COQ5, DNAJB9, MRAP, RMDN3, LPL, DNAJC15, LAMA4, SULT1A1
      22                                                                                                                                                                                                                                         IGF2, RELN, SCG3, IGFBP3, SPON1, EPHB2, SPARCL1, CA2, DOCK2, CFHR2, ITGB2, LCP1, CXCL10, IL1RL2, BTC, NGF, FLT4, ERO1A, SEMA3B, INHBA, ST6GAL1, ID2, GPNMB, IRF8, WNT7A, SATB1
      23                                                                                                                                                                                                                                                                 IL3RA, HMOX1, IL17RA, SOCS3, IL1R1, IL1R2, IL17RB, CNTFR, CRLF2, IFNGR2, CD44, IL18R1, TLR2, TNFRSF1B, CSF2RA, FAS, IRF1, IL6ST, CCL7, CSF3R, TNFRSF1A
      24                                                                                                                                                                                                    ANG, PLA2G7, PRSS3, SERPINB2, PPP4C, PRCP, CTSC, USP14, CD55, MMP8, PREP, ACTN2, PLAT, S100A9, CDA, TIMP1, CASP10, CD59, CPM, CTSD, USP15, MMP13, KYNU, ANXA5, LTF, GZMB, CFH, KCNIP3, IRF1, CTSB, ERAP2, OLR1, C1R
      25                                                                                                                                                                                                                                           ATP1B1, HSPA9, POR, TIMM13, ATP5PF, ACAT1, CYCS, PHYH, COX5B, PHB2, PDHB, COX6C, ATP5PO, GRPEL1, FXN, HTRA2, ECI1, MPC1, CS, ATP5MF, ACADVL, POLR2F, MDH2, ETFA, IDH1, COX5A
      26                                                                                                                                                                                                                                                                                                                                                                                     THY1, NRCAM, PLG, HEY1, NRP2, LDB1
      27                                                                                                                                                             IL12B, GRID2, ITIH3, BRDT, LGALS7, PCDHB1, TGM1, P2RX6, GDNF, PDCD1, DCC, EDN1, MAGIX, HNF1A, NTF3, TNNI3, SERPINA10, FGF16, ZBTB16, TGFB2, PAX4, SLC30A3, BMPR1B, LYPD3, IDUA, KCNN1, EDAR, BTG2, TG, CD80, ARPP21, AKR1B10, PTPRJ, WNT16, SLC5A5, IGFBP2
      28                                                                                                                                                                              BLVRB, CLIC3, DCXR, SORD, ALDH3B1, ZFP36, SFN, S100A9, PDLIM3, PDCD4, SULT2B1, IL17RB, UGDH, SERPINA5, CD44, COX6C, LTF, SLC27A2, KRT19, IL6ST, GLA, SERPINA3, PGR, PRKAR2B, MICB, MYB, FKBP4, CXCL14, BAG1, FGFR3, METTL3, PLAAT3, FKBP5
      29                                                                                                                                                                                                                                                                                                                                                                                        THY1, GHRL, MDGA1, NTNG1, PCSK9
      30                                                                                                                                                                                                                                                                                                                                   MB, APOD, TNNI2, TNNT2, CKM, ACHE, COL1A1, NQO1, MYOM2, LSP1, ACTN2, AK1, MYL3, ENO3
      31                                                                                                                                                                                                           ANG, PEF1, SERPINB2, MEP1A, MMP3, MMP9, APOA1, BMP1, MMP8, PREP, PLAT, KLK8, TIMP1, COMP, F9, SPARC, MST1, FURIN, CFH, CTSB, ANXA1, OLR1, C1R, MMP7, CFD, CTSH, HNF4A, DUSP6, LGMN, F8, ADAM9, SERPINE1, VWF
      32                                                                                                                                                                                                                                                                                               H2BC12, PML, HNRNPD, CUL1, HMGA1, EGF, EWSR1, CDKN1B, RPA2, ATF5, CCNA2, JPT1, CBX1, EZH2, CUL3, TOP1, TMPO, TGFB1, ODC1
      33                                                                                                                                                                                            MEP1A, SELE, IL10, CD55, NFKBIA, IL1R1, MARCO, TIMP1, CD82, NDP, CD69, CCL17, RIPK2, IFNGR2, SLAMF1, FZD5, RGS1, STAB1, TNFRSF9, IL18R1, PVR, TNFAIP6, TLR2, TNFRSF1B, IRF1, CHST2, IL10RA, CSF3, CCL7, OLR1, CSF3R, ICOSLG
      34                                                                                                                                                                                              CALR, SFN, UBE2D3, HRAS, ARF1, TBK1, EGFR, CDKN1B, UBE2N, MAP2K3, TNFRSF1A, ARHGDIA, LCK, DUSP3, PITX2, EIF4E, HSP90B1, GRK2, PLCG1, MAP2K6, PAK4, PRKCB, RIPK1, RPS6KA3, THEM4, CDKN1A, RALB, MKNK1, AKT1S1, GSK3B, GRB2
      35                                                                                                                                                                                                                                                                                                                                        CXCL10, SELL, GBP2, IRF9, PSMA3, WARS1, GMPR, IL4R, ISG15, CASP8, CXCL11, HLA-C
      36                                                                                                                                                                                                                                    THY1, CD2, IL12B, APBB1, IL13, IL12RB1, CCL22, FCGR2B, CAPG, IRF4, ITGB2, CD86, TLR3, GBP2, F2, TRAT1, IL18RAP, GZMA, WARS1, NCR1, LY75, INHBA, NCK1, TGFB2, ACVR2A, IL4R, C2, IRF8
      37                                                                                                                                                                                                                                                                                                                                                                        STC1, S100A4, TIMP1, NRP1, PRG2, APOH, SERPINA5
      38                                                                                                                                                                                                                                                                                  NOTCH2, RND3, PDGFRB, APBB2, EFEMP1, RGS4, MET, ATRN, MYC, PHF3, BMPR1A, ACVR2A, RXRA, NR1D2, DYRK1A, KIT, AGGF1, AMPH, MAGI2, INPP4B
      39                                                                                                                                                                                                                                                                                                                                                                             CRISP2, ART3, CFTR, LDHC, TSN, TSSK2, YBX2
      40                                                                                                                                                                             ALCAM, APLP1, AGER, XBP1, IL13, CA2, PENK, EEF1AKMT1, CAPG, IRF4, CXCL10, CD86, SELL, MYC, MUC1, TIAM1, IL4R, LRIG1, GLIPR2, IRF8, PIM1, PRNP, NT5E, P4HA1, CD48, BMP2, CKAP4, BMPR2, COCH, BATF3, SYT11, HOPX, TNFRSF21, RGS16, GABARAPL1
      41                                                                                                                                                  LGALS3, SOD2, CD2, LUM, ERBB3, CLU, PDGFRB, GSR, BCL2L2, F2, CTH, TGFB2, BCL2L11, CASP8, CD14, JUN, SATB1, BIK, IL18, SPTAN1, IL1A, BCL2L10, SQSTM1, BTG2, BMP2, DDIT3, BGN, EREG, NEFH, CASP7, EBP, DFFA, TNFRSF12A, FASLG, IFNB1, PTK2, BMF, BCL2L1, CD38, GSN, TNF
      42                                                                                                                                                                                                                                                                                                                                                                                                      MRTO4, NPM1, SORD
      43                                                                                                                                                                                                                                                                                                                                              NOG, LEFTY2, CDH1, BMPR1A, ID2, BMP2, SMURF1, BMPR2, SMURF2, SMAD1, BCAR3
      44                                                                                                  BPGM, PML, SOCS3, NFKBIA, ZBP1, SELP, SLAMF7, IL18BP, CD69, RIPK2, MX1, UPP1, CFH, TNFAIP6, FAS, IFIT2, IRF1, IL10RA, SAMHD1, CCL7, TXNIP, C1R, OGFR, MTHFD2, CD274, UBE2L6, BTG1, IL15RA, CD74, IL6, OASL, PNP, IL15, PSMA2, SECTM1, B2M, STAT3, IFIH1, RIPK1, LAP3, CCL5, VCAM1, SERPING1, TRIM21, NMI, PELI1, GCH1
      45                                                                                                                                                                                                                                                               ANG, GUSB, TPST1, STC1, MERTK, TGFBI, GPC1, B4GALT1, QSOX1, LDHC, GALK2, IER3, ALDOB, CHST4, HS6ST2, GLRX, EGFR, CD44, SOD1, CHST2, ALDH7A1, DCN, GAPDHS
      46                                                                                                                                                                                                                                                                         SOD2, PTPRD, NPTXR, CDKN2B, CA2, RET, CYB5R1, CCK, OLFM1, RRAD, MMP14, BCL2L11, RAB27A, SULT1A1, APOM, NR4A1, CHRNA5, SQSTM1, BTG2, HTR7, BMP2
      47                                                                                                                                                                                                                                            SERPINB2, ZFP36, SOCS3, NFKBIA, B4GALT1, GEM, CD69, IER3, RIPK2, IFNGR2, RHOB, CD44, KYNU, TNFRSF9, TNFAIP6, TLR2, IFIT2, IRF1, G0S2, TNC, IL6ST, SNN, OLR1, MAP2K3, ICOSLG
      48                                                                                                                                                                                                                                                                                                                                                                                   WNT5B, HEY1, MYC, JAG2, NOTCH1, DKK4
      49                                                                                                                                                                                                                                                                                                                        NOTCH2, DOCK2, PKD2, PLK1, VCL, ABL1, NCK1, TIAM1, BCL2L11, SPTAN1, ATG4B, LMNB1, ARHGEF2, PRC1
      50                                                                                                                                                                                                                                                                                                                                                                                                       DLL1, CUL1, FZD5
         starting_set_size final_set_size
      1                101             52
      2                201             84
      3                158            100
      4                200            119
      5                200            120
      6                150             70
      7                200            123
      8                112             41
      9                104             50
      10                40             22
      11               200            150
      12                49             31
      13               200            122
      14               200             95
      15               200            111
      16                74             41
      17               200             99
      18               200             80
      19                97             52
      20               113             57
      21               200            102
      22               201            118
      23                90             73
      24               200            146
      25               201             78
      26                36             24
      27               200             90
      28               200             96
      29                44             26
      30               200            109
      31               138            109
      32               200             79
      33               201            116
      34               105             74
      35                97             50
      36               200            140
      37                36             28
      38               144             73
      39               135             55
      40               200            109
      41               161            104
      42                58             21
      43                54             31
      44               200            117
      45               200            131
      46               158             92
      47               200            123
      48                42             22
      49               199             65
      50                32             19

---

    Code
      res$final_ranks[1:1000]
    Output
             PZP      ENPP2        NTM        LEP      ROBO2       IGF2    SLITRK1 
       14.285305   8.862328   6.520689   5.735386   5.535341   5.482352   5.426097 
       C1GALT1C1    CNTNAP2     LGALS3       SHBG  SERPINA11      ALCAM        LHB 
        5.070135   5.035036   4.997978   4.940968   4.830131   4.585340   4.542521 
        TMEM132B      CNTN1      LSAMP     GLTPD2      TIMP4      FOLR2      MDGA2 
        4.477146   4.454605   4.345116   4.200182   4.169213   4.148231   4.146763 
            PUDP        OMG    CD300LG      FETUB     ANTXR1       RELN      KCNE5 
        4.125067   4.120394   4.116513   4.113156   4.103091   4.072203   4.048407 
          ADIPOQ    COLEC12      WNT5B      LYVE1        MPZ      MGAT5      STX1A 
        3.977780   3.901979   3.881921   3.875260   3.819140   3.816185   3.802450 
            CNN1     FAM20B       THY1    C1QTNF4       A1BG       DKK2   TMEM132D 
        3.799427   3.773595   3.765527   3.763379   3.721698   3.693824   3.662841 
          SH3GL3      SEZ6L      KISS1       SOD2        CD2       GPC3      NELL2 
        3.654542   3.634272   3.630087   3.602467   3.588520   3.586675   3.581831 
           ASAH2      QPCTL       GHRL      FABP3       IAPP      OPCML      APLP1 
        3.579692   3.488765   3.442950   3.424816   3.407223   3.399488   3.370092 
             GRN       DSG2        CGA    B4GALT2       SCG3      NRCAM    SPATA24 
        3.350287   3.343559   3.343363   3.341171   3.339397   3.330785   3.312794 
           KLKB1   SERPINA7      NLGN1      MGAT2       ASPN        NMB     IGFBP3 
        3.311804   3.299856   3.298713   3.271332   3.234481   3.223505   3.212229 
            NPPB      C1QL3     CLSTN1     ATP5PB      MATN2        A2M      PILRA 
        3.181673   3.171209   3.162768   3.149317   3.141257   3.138416   3.137054 
           FABP4       MYL4    ANGPTL3      MDGA1      CILP2      BPNT2      PTPRD 
        3.131854   3.128236   3.119318   3.100271   3.099107   3.096204   3.082540 
             AGT    SLITRK5      C1QL2   C1orf162       CTSF   HEPACAM2     PCDHA4 
        3.060772   3.051341   3.049408   3.042573   3.041989   3.040368   3.031628 
      ST6GALNAC6       IL25      HABP2       VBP1     IFNA14       GNLY     ADGRB1 
        3.018906   2.997008   2.996298   2.993642   2.991857   2.989125   2.983053 
             PTN      FOXJ2    SLITRK6      KLRG2   PLA2G12B      LAIR2      CSAG2 
        2.982703   2.961016   2.946016   2.933928   2.930888   2.925945   2.922855 
         B4GALT6    C1QTNF1        OMD    TRAPPC5     CHRDL1    ANGPTL1      SPON1 
        2.916665   2.915867   2.894418   2.874215   2.870030   2.855427   2.854398 
          THSD7A        LUM        NOG     LEFTY2       CSTB       PROC      RAB26 
        2.854232   2.851378   2.844003   2.839690   2.836831   2.820480   2.820239 
          HS6ST3      SFRP4         C6        OAF      CSF1R  TNFRSF11B     LRRTM2 
        2.812797   2.808698   2.807774   2.805776   2.799469   2.788051   2.782625 
              F7      C4BPA      EPHB2      CADM2    GDAP1L1       CAST       ORC6 
        2.780071   2.775069   2.768852   2.767569   2.767150   2.764205   2.754924 
           RXFP1     FCGR2A     CRTAC1      NETO1      IL12B        UST       LSM1 
        2.738638   2.733665   2.706104   2.703544   2.701192   2.701090   2.692771 
          ZNF774      FCRL4      GRID2      TPPP2     ZNF566      PIANP      GSKIP 
        2.683923   2.683076   2.676150   2.664511   2.644861   2.643543   2.633873 
            KNG1       MYL7        CA7       BRD2      YARS2      UNC5D     SORCS1 
        2.625715   2.623557   2.621373   2.617298   2.614304   2.595381   2.593084 
            LRP8       GAS6    NAALAD2        GH2      THBS4     TMEFF1    DNAJC12 
        2.591425   2.584148   2.575922   2.563927   2.561329   2.554217   2.549145 
           DGCR6       CANX      NPTXR     LY6G6D     DNALI1     IGLON5     B3GAT1 
        2.548702   2.544511   2.541785   2.532155   2.528702   2.523679   2.521343 
         SPARCL1      MFAP4       PSG8     NOTCH2     IZUMO4      RSRP1       AGER 
        2.518862   2.518697   2.507835   2.504756   2.503449   2.499253   2.487479 
             TEF     DCTPP1       TFF3       XBP1     MAN1A2    SPATA22    S100A13 
        2.483562   2.477803   2.472687   2.468971   2.461569   2.461019   2.452950 
            PGK2      NXPH1       VCAN       NCAN       MFN1       GDF2     SLFNL1 
        2.451861   2.448173   2.446552   2.445039   2.442918   2.440365   2.439634 
         COX7A2L       CA10    SIGLEC6     SCARF2      DEPP1      CPLX1     CDKN2B 
        2.437398   2.436925   2.436438   2.434808   2.433896   2.431120   2.429321 
         ST6GAL2      KLRC1       USB1       BCL6       WWP2       MDM2        OGN 
        2.421678   2.412598   2.407627   2.407186   2.405996   2.403979   2.403220 
         ANGPTL8      NPTX1      IGLL1     IL17RE    TMEM38B       CDH1     RNF114 
        2.394512   2.392367   2.387111   2.382943   2.379199   2.368841   2.366270 
           SMYD2     LRRC15        PLG    SPATC1L      APBB1      SVEP1      ERBB3 
        2.363966   2.355264   2.352123   2.350746   2.346947   2.345193   2.339676 
           NRXN3      CHST9        CLU     HAPLN1      TBX22       MUL1       IL13 
        2.336097   2.326736   2.324165   2.323452   2.314771   2.312777   2.307899 
           ITIH3    IL12RB1        C1D        SST   C11orf49        CA2        IDS 
        2.304755   2.303107   2.299490   2.297161   2.297119   2.289490   2.288380 
          FCGR3B    CCDC126     CHST10       MREG     LRRC20      CNTN4      CLUL1 
        2.287990   2.282094   2.281336   2.281173   2.274783   2.274298   2.270427 
           S100Z       CLMP     NLGN4Y     PLXNA1      GLRX5       TLL1      SENP1 
        2.269678   2.259934   2.258965   2.258771   2.250263   2.247794   2.246684 
            PTH2        HRG        F11     PXYLP1    B3GALT2      EFNB2        DLD 
        2.241807   2.231276   2.231199   2.223169   2.207816   2.207327   2.203333 
            BRDT     UBE2J1    ADAMTS5      CBLN2    ADCYAP1       ECE1     GAGE2A 
        2.202955   2.202931   2.200561   2.197739   2.195825   2.195676   2.195338 
           EHMT2      SUSD3      TM2D1      CBLN1    FAM241B     PGRMC2     MPIG6B 
        2.193978   2.193545   2.191363   2.187880   2.180901   2.179632   2.174861 
         SULT1E1    NR2C2AP      SETD2       BMP3     HSPA1A      KLK14       IGHM 
        2.164323   2.162859   2.160223   2.159311   2.155232   2.154602   2.151964 
            CGB3    DEFB112      FCRL1      IFNA8      CRYAA      ATAD2   DEFB108B 
        2.151515   2.143189   2.142712   2.141936   2.141094   2.135267   2.135059 
           PHF11      NTNG1     CXCL17      SEMG1      SUSD1      EMC10      TREM2 
        2.132037   2.131098   2.127399   2.127141   2.124916   2.117433   2.114625 
            PENK   ADAMTSL1      GREM1      CCL22      IL36A       NDNF     FCGR2B 
        2.110253   2.105200   2.104020   2.103529   2.097636   2.093995   2.091453 
            RCN3     NIPAL4     BPIFB1    SPAG11B  EEF1AKMT1       RND3     PDGFRB 
        2.084067   2.083864   2.077421   2.077204   2.073783   2.071171   2.069777 
           APBB2       CAPG    PPP1R1A     DIXDC1      MXRA8    PLA2G2D      C1QL1 
        2.067172   2.066804   2.058576   2.058473   2.052897   2.047099   2.042815 
           CCL14         C3      PCSK1     EPHA10     HS6ST1     ZNF593      VAMP2 
        2.042124   2.042052   2.041814   2.041013   2.040950   2.040130   2.038244 
            SPRN    GABARAP      IFNA4        MAG      TCTN2    ZKSCAN7     CLEC6A 
        2.038154   2.036252   2.036083   2.035624   2.034284   2.031261   2.028585 
            SUN5        SCT      THBS2       CDSN     LGALS7      CHST3    DEFB125 
        2.028187   2.027529   2.025091   2.023839   2.021125   2.019094   2.015650 
           PSG11      VAMP3      SPIN3      DOCK2     PCDHB1        RET     HAPLN4 
        2.014509   2.012826   2.009501   2.008258   2.005848   2.005673   2.005579 
           EPHA6      CFHR2      ABCD4       PHEX      CD109       HEY1      NRIP3 
        2.005365   2.004254   2.002764   2.002367   2.001583   1.997281   1.996382 
        HS3ST3B1      PCSK9        CD5      CNTN6       FUT2      HHLA2      CD209 
        1.995520   1.993931   1.993825   1.987462   1.984027   1.982144   1.978941 
            DLG3     LRRTM4       TGM1       PKD2     BCAP29      ERP27       ARSK 
        1.971405   1.971352   1.969538   1.968500   1.962906   1.960868   1.959376 
           ITIH2       NRP2     POLR2J      ROBO1      P2RX6      NR1D1      EFNB1 
        1.958676   1.957840   1.946121   1.940018   1.939376   1.939341   1.937485 
            FGL1     IFNA10     EFEMP1    HEPACAM  TMPRSS11D      GOLM2    IL22RA2 
        1.937267   1.928237   1.927637   1.927506   1.923606   1.922239   1.921418 
          ANTXR2      DDX25     ZNF75D      CADM3      ULBP2     MRPL34   NUDT16L1 
        1.919755   1.919588   1.918123   1.914003   1.913246   1.907517   1.905009 
            KERA     NOTCH3       DTNA       SOD3       IRF4    SLC16A3      RNGTT 
        1.904735   1.902846   1.900402   1.898739   1.895394   1.891706   1.888433 
           TMEM9      ITGB2      APOA5       ASTL      OTUB2       PMP2      PTPRH 
        1.888239   1.887410   1.883026   1.881698   1.880779   1.877040   1.875769 
             LTA       GDNF      CERS5     IL11RA      PDCD1     ITGA11       FIG4 
        1.875167   1.874035   1.872956   1.871966   1.871558   1.870519   1.869911 
           NDRG3       CD33       EMC4     IL1RAP         GC    FAM172A     DYNLT3 
        1.866109   1.865417   1.858491   1.858118   1.852993   1.851778   1.850130 
            HCN1       LCP1     ANGPT2      GREM2     CXCL10       PLK1      GTF2B 
        1.849367   1.847169   1.845232   1.844544   1.843534   1.843207   1.842113 
            CD86      ARL4D       IGF1       SELL     FRRS1L     BPIFA1      KLRF1 
        1.840312   1.838732   1.834816   1.833904   1.832893   1.830452   1.829532 
            CDNF     HS3ST4       HBE1     ANKRD1      GORAB    SLITRK2     SPRED1 
        1.828436   1.822897   1.822454   1.821694   1.821264   1.818112   1.817439 
            RGS4      IFNA5      OSCAR      ITIH4       LDB1       TLR3      MAT2B 
        1.814727   1.809987   1.809846   1.809748   1.809312   1.804962   1.803027 
          CCDC51       CDY1     CYB5R1     TEPSIN  TMPRSS11A      STMN1     IGFBP1 
        1.800622   1.800603   1.797492   1.791227   1.790050   1.789012   1.788773 
          APCDD1      GIPC1        DCC      FBLN1       ARL9   TBC1D22B       WWC1 
        1.785743   1.784717   1.784655   1.784473   1.783829   1.783230   1.782879 
          BPIFA2       EDN1     TREML2      MAGIX        CCK      HNF1A       TNXB 
        1.782461   1.782268   1.778537   1.775711   1.774601   1.770449   1.769321 
          SEMA6D       C1RL     PLXDC2   SERPIND1       GDF3      TAFA5        NMU 
        1.768298   1.765796   1.763359   1.761371   1.760952   1.756891   1.755839 
            NPPA     B3GNT6     PTP4A3     RBFOX1       NTF3      RBM41      DOCK9 
        1.754566   1.753922   1.753050   1.750983   1.750661   1.749520   1.747304 
           LDOC1    MICALL2       TPPP       JPH3     COMMD5       RFX5       NRG1 
        1.743982   1.743018   1.742399   1.740840   1.739947   1.739372   1.739268 
         TNFSF12       SNX7      CSMD2       DLK1       GBP2       COQ5       SCGN 
        1.738220   1.733128   1.731137   1.730411   1.727947   1.727705   1.725767 
           FBLN7     LGALS9    EFCAB14      ILKAP      EPHB6        GSR      IL17A 
        1.725195   1.724620   1.724216   1.723811   1.722795   1.721002   1.720950 
          EGFLAM   CCNB1IP1        VCL     CT45A3     MANSC4     IL1RL2      CPLX2 
        1.719954   1.716049   1.715395   1.710115   1.709484   1.709304   1.706248 
           OPRPN      APOC1       IL19       CRNN       CDH4     DUSP18      VOPP1 
        1.705357   1.705272   1.704679   1.703021   1.701241   1.698952   1.694761 
          BCL2L2         F2      FZD10     CCDC25        BTC     DUSP16      CPXM1 
        1.693460   1.690914   1.689553   1.686515   1.686030   1.684158   1.683122 
             NGF       PSD2       POMC      FCRL5       FLT4        MET    FAM174A 
        1.680894   1.679089   1.679005   1.677492   1.676680   1.676417   1.676376 
            QPRT       SND1     DPYSL5     RNF146       FUT7   EIF4EBP1       MENT 
        1.670370   1.667837   1.666792   1.661998   1.658498   1.657744   1.657604 
           DPP10      TNNI3      UBE2T      CAMK1    TCP11L1       HMX2      SLIT2 
        1.657294   1.656721   1.656674   1.656375   1.655859   1.655521   1.654132 
         S100A11      DYRK2    DEFB119     SEMA4B    TCP10L3      UBE2C       VWC2 
        1.651932   1.646845   1.646792   1.646549   1.646448   1.643300   1.642961 
            AOC3      TRAT1       MTX2        GCG      ABCC6       ATRN       MNX1 
        1.642663   1.640332   1.634675   1.630250   1.628470   1.626402   1.625294 
          METRNL     RAB27B    IL12RB2    C1QTNF5     ADAM15      NLGN2    CRACR2A 
        1.624169   1.622504   1.620783   1.619105   1.618833   1.613872   1.609984 
           PSMA4      NMES1      RBBP6      RBBP5       WEE2     SMIM13    FAM234B 
        1.607957   1.607416   1.607256   1.606748   1.606297   1.604181   1.600468 
          GNPDA1       NRG2      FCRL6       NKD2     DNAJB9     CLEC4C       LRP1 
        1.598459   1.597367   1.595800   1.594352   1.593999   1.593868   1.592431 
           IL1RN      CRPPA      CCPG1    TMEM237      CRIP2      PCDH9        DPT 
        1.591125   1.590071   1.589225   1.588387   1.588103   1.588019   1.584924 
         IL18RAP      EMID1   TMEM167A       TRIL   SH3PXD2B       UBL3     AKR1C1 
        1.581977   1.581896   1.581853   1.578092   1.578077   1.577675   1.576439 
            BEX4  SERPINA10    GALNT13      SAR1A     ADAM30       PFN4      PRPF6 
        1.572487   1.571452   1.570311   1.569657   1.565077   1.564462   1.564272 
            FGF8    COL25A1      VAMP8      ALDOC       NRDC      NEGR1        MYC 
        1.563268   1.563152   1.562720   1.562575   1.560642   1.560460   1.557565 
            PHF3       GLCE      GDPD1       ARTN       TPM3      ETHE1      TMUB2 
        1.555831   1.555104   1.554212   1.551821   1.551212   1.550796   1.549396 
            GZMA   IL1RAPL2        INA       PALM      OLFM1      ZNRF3     PWWP2B 
        1.548117   1.544772   1.543113   1.542469   1.541323   1.539998   1.537988 
           IGHG2       DRGX      OOSP2       ADA2       CPN2      LRFN4  C10orf105 
        1.536313   1.536052   1.535219   1.533795   1.533739   1.533058   1.531922 
           FGF16       NREP       MRAP   MAP1LC3A      FNDC8     REPIN1       POLI 
        1.531527   1.530394   1.529319   1.528409   1.526768   1.524998   1.524396 
           LRIT3       OSMR      BNIP3      KCTD6      IFNA7       APTX      FDCSP 
        1.522445   1.521566   1.520687   1.519919   1.519385   1.519109   1.518982 
            CD27        HRC       EID3      MMGT1    B3GALT5   SERPINI2       CHD7 
        1.518621   1.517131   1.517031   1.516873   1.515675   1.514923   1.514652 
           MMP12        PTH       CD34   KIAA1143      BMP15       GBA3    PCOLCE2 
        1.514623   1.512254   1.511983   1.511039   1.507013   1.506509   1.505639 
          SEC61B       CHKB       PXDN       DPH5       RRAD      MMP14    SLCO5A1 
        1.504383   1.502202   1.502129   1.501997   1.501780   1.500837   1.500760 
            FGF6       ANK2      DHODH       UPRT      POSTN        CTH    WFDC10A 
        1.500006   1.499240   1.496262   1.496199   1.494712   1.493165   1.493124 
            UGP2        RP2      OBP2A    FAM171B       BRD1      ERO1A       IRF9 
        1.492425   1.491951   1.491924   1.490827   1.489551   1.489338   1.486777 
           PSMA3     SEMA3B    HSD17B7       ABL1       GBP5       MPDZ      RMDN3 
        1.486741   1.482994   1.481079   1.480525   1.480506   1.478908   1.478462 
            CA11     ZNF483     ZNF410      UBE2M      SCN2B        GHR     IGFALS 
        1.478261   1.478198   1.477630   1.475872   1.475652   1.474388   1.469903 
            MSR1     ZBTB16   LRRC37A2       CREM    ZFAND2B      WARS1      OSTM1 
        1.468701   1.467400   1.466734   1.465191   1.462795   1.462221   1.461487 
         PPP1R27      MUTYH      BTNL9      SORT1    TXNDC11      FGFR1     LILRB1 
        1.461332   1.460865   1.460124   1.460113   1.458255   1.458034   1.456325 
           APOL1       GALP        UBD    FAM102B      EGLN3        CA9      ENOX2 
        1.455926   1.454869   1.454750   1.454631   1.453850   1.451771   1.451053 
            RBL2      STX10       SSX4      FREM2      RAB2B      CALB2   SERPINA9 
        1.448004   1.447310   1.447277   1.446650   1.445111   1.444501   1.443806 
         DEFB129  TNFRSF10B    C5orf63     PGRMC1      VSIG8      MMP20       PON2 
        1.443285   1.441392   1.441062   1.440301   1.439317   1.439098   1.437781 
            ILF2      FGF23     TXNDC5     CARD19       EXOG        CCS     CRISP3 
        1.436839   1.436739   1.436406   1.435171   1.434246   1.432225   1.432079 
          ZNF580    TMEM87B      LRRC4     SPINK7    PLA2G2A     LRRC32      FANCL 
        1.428159   1.426726   1.425380   1.425055   1.424666   1.423327   1.423261 
            FUT3     IFNLR1       NCR1      FAM3B     HNRNPR      RNF34       LY75 
        1.423123   1.422346   1.420456   1.420384   1.419887   1.419798   1.419572 
            GMPR      FAM3D     MGAT4C     TOPBP1      INHBA      KCTD3       UCN3 
        1.417729   1.415012   1.412105   1.411971   1.411658   1.409796   1.408489 
          RABEPK    B3GALT1       CALY        GIP     RNF148      ARMC8     LILRA2 
        1.408120   1.407346   1.404434   1.403276   1.401178   1.400348   1.400047 
           STMN3       MUC1      LRFN3      DCP1B        CD9     CREBL2    ST3GAL1 
        1.399565   1.398395   1.398291   1.397848   1.397309   1.396633   1.396129 
           NRSN1     CYB5R3     BMPR1A      GFRAL      PSMG4      KCNE3        LPL 
        1.395409   1.394739   1.393144   1.392346   1.391278   1.391066   1.388924 
           TAFA4       RBP7       TLR4     GALNT2     POLR3F     ATP1B2       APLN 
        1.388191   1.387352   1.385941   1.385851   1.385844   1.383377   1.381990 
           BGLAP      TPPP3    COL18A1      NLRP4      FSTL5      SYT13       DDI1 
        1.380448   1.380416   1.376813   1.376237   1.376131   1.375261   1.372691 
            GRM4       DAG1     ATP1B4       SAA4    ST6GAL1      UTS2R        MOG 
        1.372566   1.372434   1.372312   1.372038   1.371468   1.370968   1.369789 
            CPN1        SF1       NCK1   HSD17B11       CAV3      ENPEP       WIF1 
        1.368695   1.367250   1.366817   1.364177   1.360203   1.360014   1.359587 
          SLC3A2       RLN3        HPD       SYT7      ASIC4       CTF1     DIRAS1 
        1.357133   1.357102   1.356578   1.356238   1.355925   1.355859   1.355576 
           TNIP1      TTC9B      CCNB1       JAG2      TRAF4       SNX1       MDM4 
        1.352841   1.351853   1.351680   1.350873   1.350777   1.350601   1.350175 
          TIMM21        DES       FAIM      TIAM1      LRIT2    BLOC1S3    HIKESHI 
        1.349876   1.349579   1.349496   1.347974   1.347311   1.347066   1.346894 
         CREB3L1     SLAMF8     CXCL16      TGFB2      SGF29       ISCU       PSG6 
        1.346880   1.346374   1.345417   1.344174   1.343997   1.343995   1.338472 
            USP2      NTRK2       TNS4       CD70      TENM3    BCL2L11    TNFSF15 
        1.338331   1.337980   1.337956   1.334391   1.328970   1.327565   1.326541 
           PSME3      PATE1     FGFBP1        CPQ       NGFR      TIMD4      CNDP1 
        1.325873   1.325290   1.323075   1.321523   1.321353   1.319685   1.319188 
            ADM2      MXRA7    SERTAD3       GRB7       KLK5       PSG1       WBP2 
        1.318739   1.317571   1.317498   1.317300   1.316531   1.316300   1.315680 
             ID2    DNAJC15     GPR101     ACVR2A   DEFB106A    KIR2DL1      PTHLH 
        1.314981   1.312473   1.312393   1.312159   1.311518   1.311137   1.310167 
             OAT       IL4R      CLCA2     DYNLT1     IGFBP4      LAMC2      CD163 
        1.309443   1.308591   1.306399   1.305389   1.305097   1.304902   1.303997 
           ISG15       PAX4    HNRNPH1      SIRT3  MACROH2A1     TCEAL7      MED11 
        1.301012   1.300317   1.291967   1.291108   1.290194   1.290079   1.290051 
           NRXN1       PES1    CNEP1R1      FLRT2     RAB27A      LRIG1     COQ10A 
        1.289910   1.289244   1.288458   1.288119   1.287936   1.287836   1.287343 
           LAMA4     MGAT4B       EURL       PSG9    LYPLAL1    PPP1R42      DHX38 
        1.287307   1.285463   1.284216   1.283265   1.282667   1.282017   1.281939 
          DNAJC4        GK5       LAD1     CD99L2     COL9A1      SNTA1       IBSP 
        1.281450   1.279330   1.278932   1.278579   1.278189   1.277657   1.276953 
            PTK7 ST6GALNAC5     LARGE1       TIE1       RXRA     GLIPR2      TOR4A 
        1.275489   1.274455   1.273927   1.272539   1.272302   1.270848   1.269534 
              C2      TFPI2      BHMT2      MYSM1      NLGN3     LIN28B     PLXDC1 
        1.267774   1.267635   1.267337   1.265289   1.264262   1.263551   1.262601 
             GBA    TXNDC12     GXYLT1   SEPTIN10     CYB5D2      RPAIN      UBE4A 
        1.262336   1.262093   1.261821   1.260902   1.260866   1.260398   1.259641 
           MZT2A     LYSMD3      TPSG1     UBXN2B       UMOD      MTMR1       ADH7 
        1.257686   1.256747   1.256522   1.256127   1.254198   1.252753   1.251972 
           CASP8    SPANXN4     LGALS2        CR2     GPR135       CD14      NPDC1 
        1.251867   1.251581   1.250317   1.249835   1.249237   1.248917   1.247307 
         JAKMIP3      RCAN1       PAK5     PCDH10      TRA2B      MYDGF      SLIT3 
        1.245552   1.243790   1.243426   1.240168   1.240057   1.239198   1.237207 
         SULT1A1    COL11A2    POMGNT2       NARF     ABHD12      POTEM      SCFD1 
        1.234871   1.233325   1.232606   1.232053   1.230847   1.230279   1.228195 
            BMT2      GPNMB       PGM5      KLRB1      ARMC5       IRF8     NT5C3B 
        1.228164   1.227736   1.227349   1.227181   1.226043   1.225391   1.225153 
           VSNL1    ADAMTS1       TCAP       AACS      NR1D2     CLEC4D   ARHGAP22 
        1.224216   1.224208   1.222934   1.220815   1.218940   1.218413   1.217776 
           PDGFD    SH3GLB1        LY9      PRRG1     B3GNT2      MCCD1      BACH1 
        1.216402   1.215354   1.214698   1.214384   1.214226   1.212412   1.212069 
            BAG3     HDGFL1       P3H1     ARRDC5      TECTB       SUOX    FAM204A 
        1.211390   1.211262   1.211087   1.210750   1.210260   1.210056   1.209439 
             AMH    SLC30A3       MZF1    MACROD2     KLHL40      FGFR4     ADAM22 
        1.208100   1.207509   1.207303   1.206809   1.206124   1.204479   1.203406 
           KRT20      DDIT4  SERPINA12      WNT7A        JUN      SATB1       SGK3 
        1.203052   1.202299   1.200432   1.199649   1.198340   1.197592   1.196330 
            POLH      CRIM1      PDE1B      ESRP1        BIK      VAMP4      HGFAC 
        1.193615   1.192993   1.192321   1.191761   1.191743   1.190044   1.189672 
         CCDC134       APOM      ECRG4      CRYGS       TSR2      ICAM2     ZFAND3 
        1.189214   1.186697   1.184849   1.184480   1.183993   1.182258   1.181795 
           POLD4       SVIP      KRT72    BLOC1S5     ACVR2B    ALDH5A1       CLK2 
        1.180560   1.179157   1.178697   1.177180   1.176152   1.175829   1.175159 
           IFNA6       DLX3    KIRREL1      TEAD4     EXOSC3       PIM1 
        1.174775   1.174052   1.173488   1.173302   1.171418   1.171339 

# `somaPrGSEA(verbose = TRUE)` emits expected messages

    Code
      t <- somaPrGSEA(ranks = ranks_noDimer, verbose = TRUE)
    Output
      -- Resolving Non-Unique Mapping ------------------------------------------------
    Message
      x Duplicate value(s) found in the ranking vector names.
      > Filtering will be performed to retain a single value.
      > There are 772 feature(s) with multiple elements in the ranked input vector.
      v A total of 877 lower-ranked element(s) were dropped.
    Output
      -- Performing Enrichment Analysis ----------------------------------------------
    Message
      > Running pre-ranked GSEA...
      i If `split_ids = TRUE`, any heterodimers present were split into individual elements.
      i This may introduce ties. See `?somaPrGSEA()`.
      v Done!

