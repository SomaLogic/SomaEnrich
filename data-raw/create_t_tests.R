library(dplyr)
library(purrr)
library(SomaDataIO)
library(tidyr)

# Code below taken from SomaDataIO "Two-Group Comparison" workflow and
# modified to use 11K ADAT as input.

# prepare data set for analysis using `preProcessAdat()`
cleanData <- example_data_11k |>
    preProcessAdat(
        filter.features = TRUE,            # rm non-human protein features
        filter.controls = TRUE,            # rm control samples
        filter.qc       = TRUE,            # rm non-passing qc samples
        log.10          = TRUE,            # log10 transform
        center.scale    = FALSE            # do not center/scale analytes
    ) |>
    drop_na(Sex)                           # rm NAs if present

table(cleanData$Sex, useNA = "always")

t_tests <- getAnalyteInfo(cleanData)

any(is.na(t_tests$EntrezGeneSymbol))

t_tests <- t_tests |>
    filter(!EntrezGeneSymbol %in% c("None", "")) |>        # Addl step for this R pkg, not in original workflow
    select(AptName, SeqId, Target, EntrezGeneSymbol, UniProt) |>
    mutate(
        formula = map(AptName, ~ as.formula(paste(.x, "~ Sex"))), # create formula
        t_test  = map(formula, ~ stats::t.test(.x, data = cleanData)),  # fit t-tests
        t_stat  = map_dbl(t_test, "statistic"),            # pull out t-statistic
        p.value = map_dbl(t_test, "p.value"),              # pull out p-values
        fdr     = p.adjust(p.value, method = "BH")         # FDR for multiple testing
    )

all_seqs <- getAnalytes(cleanData)
all_seqs <- intersect(all_seqs, t_tests$AptName)

# Calculate FC for each analyte (columns) using M as baseline.
# Match by AptName so row order of t_tests doesn't matter.
fc <- apply(cleanData[, all_seqs], 2, function(x) {
    median_case <- median(x[cleanData$Sex == "F"])
    median_ctrl <- median(x[cleanData$Sex == "M"])
    median_case / median_ctrl
})

t_tests$fc <- fc[t_tests$AptName]
t_tests$log2_fc <- log2(t_tests$fc)
# t_tests$neglog10p <- -log10(t_tests$p.value)*sign(t_tests$log2_fc)

t_tests <- t_tests |>
    arrange(desc(t_stat), fdr)

save(t_tests, file = "data/t_tests.rda", compress = "xz")
