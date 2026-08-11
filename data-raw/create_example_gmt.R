# Download example GMT from MSigDb.
# This is the Reactome subset of the canonical pathways subcollection (in C2)
gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c2.cp.reactome.v2025.1.Hs.symbols.gmt"
download_loc <- "inst/extdata/msigdb_c2_reactome.gmt"

ex_gmt <- download.file(gmt_url, destfile = download_loc)

# Check that procedure above worked as expected
test <- fgsea::gmtPathways(gmt.file = download_loc)

head(test)
