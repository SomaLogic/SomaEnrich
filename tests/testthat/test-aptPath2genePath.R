# Setup ------

# This step is important, not all checks will pass without it
# TODO: do the tests need to be modified to not require this step, to be 
# more robust?
anno <- getAnalyteInfo(example_data_11k) |>
  tidyr::separate_rows(EntrezGeneSymbol, sep = "\\,|\\s+|\\|")

# Differs from example_pathway because all values are found in extended annots
pthwy_complete <- withr::with_seed(123, {
  anno |>
    dplyr::filter(EntrezGeneSymbol != "") |>
    dplyr::pull(AptName) |>
    sample(50, replace = FALSE)
})

all_pthwy <- list(complete = pthwy_complete,
                  missing = ex_apt_pathway)


# Testing -----
test_that("`aptPath2genePath()` produces expected output with default args", {
  out <- aptPath2genePath(path_list = all_pthwy, col_meta_df = anno)
  
  expect_length(out, 2L)
  expect_type(out, "list")
  expect_false(any(duplicated(names(out))))
  
  # Make sure each item is unique
  expect_false(identical(out$complete, out$missing))
  expect_length(out$complete, 51L)
  expect_length(out$missing, 58L) # TODO: Test to make sure expected genes are missing
  
  # Requires `separate_rows(anno)` call from setup
  expect_true(all(out[[1]] %in% anno$EntrezGeneSymbol)) 
})

test_that("`aptPath2genePath()` handles Entrez IDs", {
  out <- aptPath2genePath(path_list = all_pthwy[1], 
                          col_meta_df = anno, 
                          id_type = "EntrezGeneID")
  
  # Don't have exact same data for symbols and IDs, this is expected to be 
  # different from when symbols are used
  expect_length(out, 1L)
  expect_type(out, "list") 
  expect_length(out$complete, 51L)
  expect_type(out$complete, "character") 
  
  expect_false(any(duplicated(out$complete)))
  
  # Requires `separate_rows(anno)` call from setup
  expect_true(all(out$complete %in% anno$EntrezGeneID)) 
})

test_that("`aptPath2genePath()` catches when non-AptNames are provided", {
    expect_error(
        aptPath2genePath(path_list = list(fail_path = c("12345-67", "seq.98765.43", "test")), 
                         col_meta_df = anno,
                         verbose = TRUE),
        "All SomaScan identifiers should be in 'AptName' format."
    )
})


test_that("`aptPath2genePath(verbose = TRUE)` works as expected", {
  expect_message(
    aptPath2genePath(path_list = list(msg_path = c("seq.12345.67", "seq.98765.43", "seq.2953.31")), 
                     col_meta_df = anno,
                     verbose = TRUE),
    "2 analytes in 'path_list' do not have an associated gene identifier"
  )
  
  expect_silent(
    aptPath2genePath(path_list = all_pthwy[1], 
                     col_meta_df = anno,
                     verbose = FALSE)
  )
})

test_that("`aptPath2genePath()` returns all genes when 1:many mapping is present", {
  apt <- "seq.5245.40"
  out <- aptPath2genePath(path_list = list(test = apt), col_meta_df = anno)
  
  expect_length(out$test, 3L)
  expect_true(all(c("PRKAA2", "PRKAB2", "PRKAG1") %in% out$test))
})

test_that("`aptPath2genePath()` output has no duplication when 1:1 mapping is present", {
  # All map to gene CGA
  cga_apts <- c("seq.4914.10", "seq.3032.11", "seq.2953.31")
  out <- aptPath2genePath(path_list = list(test = cga_apts), col_meta_df = anno)
  
  expect_length(out$test, 5L)
  expect_false(any(duplicated(out$test)))
})

test_that("`aptPath2genePath()` throws error when non-AptNames are used as input", {
  expect_error(
    aptPath2genePath(path_list = list(test = c("1234-56", "seq.1234.56", "seq.5678.90")), 
                     col_meta_df = anno),
    "All SomaScan identifiers should be in 'AptName' format."
  )
})

test_that("`aptPath2genePath()` works when col_meta_df has dot-separated column name", {
  # Create dummy IPP (NGS) ADAT metadata
  anno_ipp <- anno
  names(anno_ipp)[names(anno_ipp) == "EntrezGeneSymbol"] <- "Entrez.Gene.Symbol"
  
  res_ipp <- aptPath2genePath(path_list = all_pthwy[1], col_meta_df = anno_ipp)
  res_array <- aptPath2genePath(path_list = all_pthwy[1], col_meta_df = anno)
  
  # Make sure results are the same, regardless of NGS/array ADAT
  expect_equal(res_ipp, res_array)
})

test_that("`aptPath2genePath()` throws error when wrong 'id_type' used", {
  expect_error(
    aptPath2genePath(path_list = all_pthwy[1], 
                     col_meta_df = anno, 
                     id_type = "UniProt"),
    "'arg' should be one of \"EntrezGeneSymbol\", \"EntrezGeneID\""
  )
})

test_that("`aptPath2genePath()` returns empty gene vectors for pathways with no matching analytes", {
  # All analytes are valid AptNames but not in col_meta_df
  missing_path <- list(NoMatch = c("seq.12345.67", "seq.98765.43"))
  out <- aptPath2genePath(path_list = missing_path,
                           col_meta_df = anno,
                           verbose = FALSE)
  
  expect_length(out, 1L)
  expect_length(out$NoMatch, 0L)
  expect_type(out$NoMatch, "character")
})

test_that("`aptPath2genePath()` silently passes when path contains an empty AptName vector", {
  # Known limitation: all(is.AptName(unlist(list(x = character(0))))) -> TRUE
  # because all(logical(0)) is TRUE in R. This documents the current behavior.
  empty_path <- list(EmptyPath = character(0))
  
  # Should NOT error — empty unlist passes the is.AptName() check silently
  expect_no_error(
    aptPath2genePath(path_list = empty_path, col_meta_df = anno, verbose = FALSE)
  )
  out <- aptPath2genePath(path_list = empty_path, col_meta_df = anno, verbose = FALSE)
  expect_length(out$EmptyPath, 0L)
})

test_that("`aptPath2genePath()` handles a path_list with multiple entries, some empty", {
  mixed_paths <- list(
    Valid   = pthwy_complete[1:5],
    Empty   = character(0),
    Missing = c("seq.12345.67")  # Valid AptName format but not in col_meta_df
  )
  out <- aptPath2genePath(path_list = mixed_paths,
                           col_meta_df = anno,
                           verbose = FALSE)
  
  expect_length(out, 3L)
  expect_gt(length(out$Valid), 0L)    # Valid apts should map to genes
  expect_length(out$Empty, 0L)        # Empty input -> empty output
  expect_length(out$Missing, 0L)      # Unknown apt -> empty output
})
