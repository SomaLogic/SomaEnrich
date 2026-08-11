# ------ Setup
adat <- SomaDataIO::example_data

df_list <- SomaDataIO::rn2col(adat, "Rowname") |>
  df2list(name_col = "SampleType", value_col = "Rowname")


# ------ Testing
test_that("`df2list()` produces expected result with default arguments", {
  expect_type(df_list, "list")
  expect_named(df_list, c("Sample", "Calibrator", "Buffer", "QC"))
  expect_identical(df_list$QC, c("258495800001_6", "258495800007_6",
                                 "258495800002_3", "258495800109_6",
                                 "258495800112_6", "258495800107_7")
  )
})

test_that("`df2list()` errors if name and value columns aren't present", {
  expect_error(
    df2list(adat[SomaDataIO::getMeta(adat)],
            name_col = "SampleType", value_col = "12345"),
    "`value_col` must be a column name"
  )
  expect_error(
    df2list(adat[SomaDataIO::getMeta(adat)],
            name_col = "12345", value_col = "SampleId"),
    "`name_col` must be a column name"
  )
})

test_that("`df2list()` works on grouped data frames", {
  grp_df <- adat |>
    dplyr::group_by(PlatePosition)
  l <- df2list(grp_df, name_col = "SampleType", value_col = "PlatePosition")
  
  expect_named(l, c("Sample", "Calibrator", "Buffer", "QC"))
  expect_equal(l$Buffer, c("G9", "E3", "A7", "E10", "D3", "B5" ))
})
