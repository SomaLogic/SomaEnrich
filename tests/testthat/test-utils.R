
test_that(".is.sorted() returns TRUE for ascending sorted vector", {
    expect_true(.is.sorted(c(1, 2, 3, 4, 5)))
})

test_that(".is.sorted() returns TRUE for descending sorted vector", {
    expect_true(.is.sorted(c(5, 4, 3, 2, 1)))
})

test_that(".is.sorted() returns FALSE for unsorted vector", {
    expect_false(.is.sorted(c(3, 1, 4, 1, 5)))
})

test_that(".is.sorted() returns TRUE or FALSE (never NA) when NAs are present", {
    na_asc  <- c(1, 2, NA, 4)   # ascending with NA
    na_desc <- c(4, NA, 2, 1)   # descending with NA
    na_unsorted <- c(3, NA, 1)  # unsorted with NA
    
    expect_true(!is.na(.is.sorted(na_asc)))
    expect_true(!is.na(.is.sorted(na_desc)))
    expect_true(!is.na(.is.sorted(na_unsorted)))
})


####### .match_col() ########

test_that(".match_col() returns column name unchanged on exact match", {
    df <- data.frame(EntrezGeneSymbol = 1, AptName = 2)
    expect_equal(.match_col(df, "EntrezGeneSymbol"), "EntrezGeneSymbol")
    expect_equal(.match_col(df, "AptName"), "AptName")
})

test_that(".match_col() matches dot-separated variant (Entrez.Gene.Symbol)", {
    df <- data.frame(Entrez.Gene.Symbol = 1, AptName = 2)
    expect_equal(.match_col(df, "EntrezGeneSymbol"), "Entrez.Gene.Symbol")
})

test_that(".match_col() matches all-lowercase variant", {
    df <- data.frame(entrezgenesymbol = 1, AptName = 2)
    expect_equal(.match_col(df, "EntrezGeneSymbol"), "entrezgenesymbol")
})

test_that(".match_col() matches EntrezGeneID dot-separated variant", {
    df <- data.frame(Entrez.Gene.ID = 1, AptName = 2)
    expect_equal(.match_col(df, "EntrezGeneID"), "Entrez.Gene.ID")
})

test_that(".match_col() errors informatively when no match found", {
    df <- data.frame(SomeOtherColumn = 1, AptName = 2)
    expect_error(
        .match_col(df, "EntrezGeneSymbol"),
        "No column matching 'EntrezGeneSymbol' found in the provided data frame."
    )
})

test_that(".match_col() errors informatively on ambiguous match", {
    # Two columns that normalize to the same string
    df <- data.frame(Entrez.Gene.Symbol = 1, Entrez_Gene_Symbol = 2, AptName = 3)
    expect_error(
        .match_col(df, "EntrezGeneSymbol"),
        "Ambiguous column match for 'EntrezGeneSymbol'"
    )
})

test_that(".match_col() error propagates correctly through public functions", {
    meta <- getAnalyteInfo(example_data_11k)
    meta_noGene <- data.frame(AptName = meta$AptName, 
                              TestCol = 1)
    expect_error(
        apt2gene("seq.1111.11", col_meta_df = meta_noGene),
        "No column matching 'EntrezGeneSymbol' found in the provided data frame."
    )
    expect_error(
        gene2apt("Gene1", col_meta_df = meta_noGene),
        "No column matching 'EntrezGeneSymbol' found in the provided data frame."
    )
})
