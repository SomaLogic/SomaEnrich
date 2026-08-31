# Select Best Value Based on Metric

When multiple values map to the same protein/gene identifier, retains a
single value to represent the gene based on a pre-defined selection
method.

## Usage

``` r
.select_best(
  ranks,
  resolve_method = c("abs", "min", "max", "rank"),
  call = rlang::caller_env()
)
```

## Arguments

- ranks:

  A named numeric vector.

- resolve_method:

  Character. Selection method for choosing a single representative value
  for a gene when non-1:1 mapping exists between analytes and genes.
  Choose one of:

  "abs"

  :   (Default) Selects the entry with the largest absolute value per
      gene. Recommended for signed statistics (e.g. t-statistic, log
      fold change) in bidirectional analyses — unlike `"rank"` or
      `"max"`, it preserves strong signal symmetrically at both the top
      and bottom of the ranking vector.

  "min"

  :   Selects the entry with the minimum numeric value, including
      negative values. Use for metrics where smaller values are better
      (e.g., p-values) or when smaller/negative signed values are
      important.

  "max"

  :   Selects the entry with the maximum numeric value. Use for metrics
      where larger values are preferred (e.g. positive one-sided tests).

  "rank"

  :   Selects the entry with the smallest positional index (i.e. first
      occurrence in the input vector). Use when input is pre-sorted by
      one or more variables, and the "best" entry is the one closest to
      position 1. This means that the selected entry depends on the
      order of the input, so `ranks` must be **pre-sorted** to place
      values of highest importance first.

## Value

List containing:

- ranks:

  Ranking vector used to select features. May differ from input ranks if
  heterodimers are present, as they will be split into individual
  components.

- id_list:

  List of features. Elements are values from the `ranks` vector.

- keep_list:

  List of features, with only one value per feature. List elements are
  the index of the retained value (based on "method") from the original
  `ranks` vector.
