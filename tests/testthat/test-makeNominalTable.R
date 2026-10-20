# Setup ------
contrasts <- data.frame(group    = c("A", "A", "B", "B"),
                        estimate = c(0.8, 0.1, -0.9, 0.6),
                        p.value  = c(0.01, 0.01, 0.001, 0.2))


# Testing -----
test_that("`makeNominalTable()` rejects reserved `group_col` names", {
  for ( nm in c("cutoff", "fc", "Count") ) {
    df <- contrasts
    df[[nm]] <- df$group
    expect_error(
      makeNominalTable(df, group_col = nm),
      "`group_col` cannot be one of the reserved output column names"
    )
  }
})

test_that("`makeNominalTable()` requires `group_col` to be a single string", {
  expect_error(
    makeNominalTable(contrasts, 
                     group_col = c("group", "estimate")),
    "`group_col` must be a single character string."
  )
})

test_that("`makeNominalTable()` still works with a valid `group_col`", {
  out <- makeNominalTable(contrasts, group_col = "group",
                          sig.levels = 0.05, fc.levels = 0.5)
  
  expect_named(out, c("group", "cutoff", "fc", "Count"))
  expect_equal(out$group, c("A", "B"))
  expect_equal(out$Count, c(1L, 1L))
})

test_that("`makeNominalTable()` counts every cutoff pair within groups", {
    out <- makeNominalTable(contrasts, group_col = "group",
                            sig.levels = c(0.01, 0.2),
                            fc.levels = c(0.1, 0.8, 0.9))
    expected <- data.frame(group  = rep(c("A", "B"), 6),
                           cutoff = rep(c(0.01, 0.2), each = 6),
                           fc     = rep(rep(c(0.1, 0.8, 0.9), each = 2), 2),
                           Count  = c(2L, 1L, 1L, 1L, 0L, 1L, 2L, 2L, 1L, 1L, 0L, 1L))
    
    expect_equal(out, expected)
})
test_that("`makeNominalTable()` counts every cutoff pair without grouping", {
    out <- makeNominalTable(contrasts, 
                            sig.levels = c(0.01, 0.2),
                            fc.levels  = c(0.1, 0.8, 0.9))
    expected <- data.frame(cutoff = rep(c(0.01, 0.2), each = 3),
                           fc     = rep(c(0.1, 0.8, 0.9), 2),
                           Count  = c(3L, 2L, 1L, 4L, 2L, 1L))
    
    expect_equal(out, expected)
})

test_that("`makeNominalTable()` returns a zero-row table for empty cutoffs", {
  out <- makeNominalTable(contrasts, sig.levels = numeric(0))
  
  expect_s3_class(out, "data.frame")
  expect_named(out, c("cutoff", "fc", "Count"))
  expect_equal(nrow(out), 0L)
})
