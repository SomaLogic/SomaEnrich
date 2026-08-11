# ------ Setup
meta <- SomaDataIO::getAnalyteInfo(example_data_11k)
gene_df <- collapseAptData(df = meta)


# ------ Testing
test_that("`collapseAptData()` works as expected with default parameters", {
    unique_genes <- meta |>
        tidyr::separate_rows(EntrezGeneSymbol, sep = "\\,|\\s+|\\|") |>
        dplyr::filter(EntrezGeneSymbol != "") |>
        dplyr::pull(EntrezGeneSymbol) |>
        unique()
    apt_order <- gene_df[gene_df$EntrezGeneSymbol == "APOE", "AptName"]
    
    expect_s3_class(gene_df, "data.frame")
    expect_length(gene_df$EntrezGeneSymbol, length(unique_genes))
    expect_equal(apt_order, "seq.2418.55|seq.2937.10|seq.2938.55|seq.5312.49")
})

test_that("`collapseAptData()` can use different delimiters for collapsed rows", {
    res_delim <- collapseAptData(meta, sep = ", ")
    ex_collapsed <- res_delim[res_delim$EntrezGeneSymbol == "KLK3", ]$AptName
    expect_true(grepl("\\, ", ex_collapsed))
})

test_that("`collapseAptData()` reorders analytes when `apt_order = TRUE`", {
    # Default order of APOE-associated SeqIds in dataset
    def_apoe <- dplyr::filter(meta, EntrezGeneSymbol == "APOE")$AptName
    
    # Creating an arbitrary, random order
    meta$Analyte_Rank <- sample(1:nrow(meta), nrow(meta), replace = FALSE)
    
    # New, modified order of APOE SeqIds
    newOrd_apoe <- dplyr::arrange(meta, Analyte_Rank) |> 
        dplyr::filter(EntrezGeneSymbol == "APOE") |>
        dplyr::pull(AptName) |>
        paste(collapse = "|")
    
    res_ord <- collapseAptData(meta, apt_order = "Analyte_Rank")
    res_apoe <- res_ord[res_ord$EntrezGeneSymbol == "APOE", "AptName"]
    
    expect_equal(res_apoe, newOrd_apoe)  
})

test_that("`collapseAptData()` generates a list-column when sep = NULL", {
    list_df <- collapseAptData(meta, sep = NULL)
    
    expect_type(list_df$AptName, "list")
    expect_equal(list_df$AptName$ABL2, c("seq.3342.76", "seq.5261.13"))
})

test_that("`collapseAptData()` works when col_meta_df has dot-separated column name", {
    meta_ipp <- meta
    names(meta_ipp)[names(meta_ipp) == "EntrezGeneSymbol"] <- "Entrez.Gene.Symbol"
    res_ipp   <- collapseAptData(meta_ipp)
    res_array <- collapseAptData(meta)
    # Output column name should reflect the actual name found in the data frame
    expect_named(res_ipp, c("Entrez.Gene.Symbol", "AptName"))
    expect_equal(res_ipp[["Entrez.Gene.Symbol"]], res_array[["EntrezGeneSymbol"]])
    expect_equal(res_ipp[["AptName"]], res_array[["AptName"]])
})

test_that("`collapseAptData()` with sep = NULL works with dot-separated column name", {
    meta_ipp <- meta
    names(meta_ipp)[names(meta_ipp) == "EntrezGeneSymbol"] <- "Entrez.Gene.Symbol"
    list_ipp   <- collapseAptData(meta_ipp, sep = NULL)
    list_array <- collapseAptData(meta, sep = NULL)
    expect_named(list_ipp, c("Entrez.Gene.Symbol", "AptName"))
    expect_equal(list_ipp[["AptName"]], list_array[["AptName"]])
})

test_that("`collapseAptData()` warns when `apt_order` column is non-numeric", {
    meta_char <- meta
    meta_char$CharRank <- sample(letters, nrow(meta_char), replace = TRUE)
    
    expect_warning(
        collapseAptData(meta_char, apt_order = "CharRank"),
        "Ordering may produce unexpected results when the `apt_order` column is not numeric."
    )
})

test_that("`collapseAptData()` handles a data frame with a single row", {
    single_row <- meta[1L, ]
    res <- collapseAptData(single_row)
    
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 1L)
    expect_named(res, c("EntrezGeneSymbol", "AptName"))
})

test_that("`collapseAptData()` handles a data frame where all genes map 1:1", {
    # Build a truly 1:1 mock (one AptName per gene, no multi-gene strings)
    meta_11 <- data.frame(
        AptName          = c("seq.1111.11", "seq.2222.22", "seq.3333.33"),
        EntrezGeneSymbol = c("GeneA",       "GeneB",       "GeneC")
    )
    res <- collapseAptData(meta_11)
    
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 3L)
    # No pipes: every AptName cell should be a plain AptName, no collapsing needed
    expect_false(any(grepl("\\|", res$AptName)))
})
