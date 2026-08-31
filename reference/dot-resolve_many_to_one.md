# Select A Single Measurement Per Feature

Select A Single Measurement Per Feature

## Usage

``` r
.resolve_many_to_one(
  ranks,
  resolve_method = c("abs", "min", "max", "rank"),
  verbose = interactive()
)
```

## Arguments

- ranks:

  Either the output of
  [`prepareRanks()`](https://somalogic.github.io/SomaEnrich/reference/prepareRanks.md),
  or a named vector of a numeric ranking statistic. This vector will be
  used to rank features, according to the value of the statistic. Any
  numeric metric may be used (t-statistic, -log10(p), fold change,
  etc.). The vector *must* be named using Entrez gene
  symbols/identifiers *or* SomaScan identifiers in R-compatible
  `AptName` format (e.g. `seq.1234.56`). If using `AptNames`, set
  `use_aptnames = TRUE`.

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

- verbose:

  Logical. Should informative progress messages be printed to the
  console?

## Value

A list containing two elements:

- results:

  The filtered input vector, with non-prioritized features removed.

- removed:

  A named vector containing elements removed from the input `ranks`
  vector.
