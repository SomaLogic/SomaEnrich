# Setup ----

# Test strings
named_pipe <- c("GeneA|GeneB" = 1.5, "GeneC" = 2.0, "GeneD GeneE GeneF" = 3.0)
unnamed_pipe <- names(named_pipe)
named_comma <- c("GeneA,GeneB" = 1.5, "GeneC" = 2.0)
named_single <- c("GeneA" = 1.0, "GeneB" = 2.0, "GeneC" = 3.0)

# "Real world" example:
# Total unique groupings of genes & AptNames in the ex differential exp results
genes_7k <- t_tests |>
    dplyr::select(AptName, EntrezGeneSymbol) |>
    tidyr::separate_rows(EntrezGeneSymbol, sep = "\\,|\\s+|\\|")

# Reference values to test multimer filtering
n_split_genes <- length(genes_7k$EntrezGeneSymbol)
n_unique_genes <- length(unique(genes_7k$EntrezGeneSymbol))



# Testing ----
test_that("`.splitIDs()` correctly splits named vector with default delimiter", {
    pipe_res <- .splitIDs(named_pipe)
    
    expect_type(pipe_res, "double")
    expect_named(pipe_res)
    expect_length(pipe_res, 6)  # 2 + 1 + 3 = 6 total, after splitting
    expect_equal(names(pipe_res), c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"))
    expect_equal(unname(pipe_res), c(1.5, 1.5, 2.0, 3.0, 3.0, 3.0))
    expect_identical(unname(pipe_res["GeneA"]), unname(pipe_res["GeneB"]))
})

test_that("`.splitIDs()` correctly splits named vector with user-specified delimiter", {
    result <- .splitIDs(named_comma, ",")
    
    expect_type(result, "double")
    expect_named(result)
    expect_length(result, 3)  # 2 + 1 = 3 total entries
    expect_equal(names(result), c("GeneA", "GeneB", "GeneC"))
    expect_equal(unname(result), c(1.5, 1.5, 2.0))
})

test_that("`.splitIDs()` correctly splits unnamed vector with default delimiter", {
    result <- .splitIDs(unnamed_pipe)
    
    expect_type(result, "character")
    expect_null(names(result))
    expect_length(result, 6L)  # 2 + 1 + 3 = 6 total entries
    expect_equal(result, c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"))
})

test_that("`.splitIDs()` handles vector with no delimiters (pass-through)", {
    result <- .splitIDs(named_single)
    
    expect_type(result, "double")
    expect_named(result)
    expect_length(result, 3L)
    expect_equal(names(result), c("GeneA", "GeneB", "GeneC"))
    expect_equal(unname(result), c(1.0, 2.0, 3.0))
})

test_that("`.splitIDs()` handles mixed delimiters in unnamed vector", {
    mixed_delim <- c("GeneA|GeneB", "GeneC,GeneD", "GeneE GeneF")
    result <- .splitIDs(mixed_delim)
    
    expect_type(result, "character")
    expect_length(result, 6)
    expect_equal(result, c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"))
})

test_that("`.splitIDs()` handles single-element input", {
    single_named <- c("GeneA" = 5.0)
    single_unnamed <- "GeneA"
    
    result_named <- .splitIDs(single_named)
    result_unnamed <- .splitIDs(single_unnamed)
    
    expect_length(result_named, 1)
    expect_equal(names(result_named), "GeneA")
    expect_equal(unname(result_named), 5.0)
    
    expect_length(result_unnamed, 1)
    expect_equal(result_unnamed, "GeneA")
})

test_that("`.splitIDs()` handles single element input with delimiter", {
    single_named <- c("GeneA|GeneB" = 5.0)
    single_unnamed <- "GeneA|GeneB"
    
    result_named <- .splitIDs(single_named)
    result_unnamed <- .splitIDs(single_unnamed)
    
    expect_length(result_named, 2)
    expect_equal(names(result_named), c("GeneA", "GeneB"))
    expect_equal(unname(result_named), c(5.0, 5.0))
    
    expect_length(result_unnamed, 2)
    expect_equal(result_unnamed, c("GeneA", "GeneB"))
})

test_that("`.splitIDs()` preserves value types", {
    # Integer values
    int_named <- c("GeneA|GeneB" = 1L, "GeneC" = 2L)
    result_int <- .splitIDs(int_named)
    expect_type(result_int, "integer")
    
    # Character values in unnamed vector
    char_input <- c("GeneA|GeneB", "GeneC")
    result_char <- .splitIDs(char_input)
    expect_type(result_char, "character")
})

test_that("`.splitIDs()` produces expected output with SomaScan example data", {
    t_test_res <- .splitIDs(t_tests$EntrezGeneSymbol)
    
    expect_length(t_test_res, n_split_genes)
    expect_length(unique(t_test_res), n_unique_genes)
})
