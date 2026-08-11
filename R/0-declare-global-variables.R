# --------------------------- #
# Declaring Global Variables:
# This is mostly for passing R CMD checks
# global variables that come from other dependant
# packages, or objects in the 'data/' directory
# Reference: https://github.com/tidyverse/magrittr/issues/29
# ---------------------------------------------------------- #
if ( getRversion() >= "2.15.1" )
  utils::globalVariables(
    c(".data",
      "category",
      "example_data_11k",
      "foldEnrichment",
      "geneList",
      "hit_color",
      "id",
      "overlap",
      "padj",
      "pctDE",
      "ranks",
      "runningScore",
      "pathway_map",
      "UniProt",
      "x")
  )
