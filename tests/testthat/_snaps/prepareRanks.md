# `prepareRanks(verbose = TRUE)` emits expected messages

    Code
      t <- prepareRanks(stats = dup_stats, features = dup_features, split_ids = FALSE,
        resolve_multimapping = TRUE, verbose = TRUE)
    Output
      -- Resolving Non-Unique Mapping ------------------------------------------------
    Message
      x Duplicate value(s) found in the ranking vector names.
      > Filtering will be performed to retain a single value.
      > There are 1 feature(s) with multiple elements in the ranked input vector.
      v A total of 1 lower-ranked element(s) were dropped.

# `prepareRanks()` warns about NAs, with the correct count

    Code
      prepareRanks(stats = na_stats, features = na_features, verbose = FALSE)
    Condition
      Warning:
      2 NA values detected in `stats`. They may be dropped if `resolve_multimapping = TRUE`.
    Output
      GeneA GeneB GeneC GeneD 
         NA     2    NA     4 

---

    Code
      prepareRanks(stats = na_stats[-1], features = na_features[-1], verbose = FALSE)
    Condition
      Warning:
      1 NA value detected in `stats`. They may be dropped if `resolve_multimapping = TRUE`.
    Output
      GeneB GeneC GeneD 
          2    NA     4 

