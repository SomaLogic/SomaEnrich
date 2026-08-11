library(msigdbr)
library(dplyr)
library(SomaDataIO)

colxn <- "H"
path_name <- "HALLMARK_ANDROGEN_RESPONSE"

ex_gene_pathway <- msigdbr::msigdbr(collection = colxn) |>
    dplyr::filter(gs_name == path_name) |>
    dplyr::pull(gene_symbol) |>
    unique()

ex_apt_pathway <- gene2apt(x = ex_gene_pathway, col_meta_df = SomaDataIO::getAnalyteInfo(example_data_11k), collapse = FALSE)
rmv <- which(!is.na(ex_apt_pathway))
ex_apt_pathway <- unique(ex_apt_pathway[!is.na(ex_apt_pathway)])
ex_gene_pathway <- ex_gene_pathway[-rmv]

# Checks
any(duplicated(ex_gene_pathway))
any(duplicated(ex_apt_pathway))

save(ex_gene_pathway, file = "data/ex_gene_pathway.rda", compress = "xz")
save(ex_apt_pathway, file = "data/ex_apt_pathway.rda", compress = "xz")
