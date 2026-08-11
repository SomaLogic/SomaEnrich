suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(dplyr)
  library(GO.db)
  library(stringr)
  library(msigdbr)
  library(org.Hs.eg.db)
})


################################
#--- Retrieve pathway data ----#
################################

###### Gene Ontology (GO) #########################

go_term_map <- AnnotationDbi::Ontology(GO.db::GOTERM)
terms2keep  <- go_term_map[go_term_map %in% c("BP", "MF")]
go_term_map <- data.frame(term   = names(go_term_map),
                          aspect = unname(go_term_map)) # Will use later

go_name_map <- AnnotationDbi::keys(GO.db, keytype = "GOID")
go_name_map <- AnnotationDbi::select(GO.db, columns = c("GOID", "TERM"),
                                     keys = go_name_map, keytype = "GOID")

# Important to suppress messages! Very slow otherwise
go_symbol_list <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(terms2keep),
                          column = "SYMBOL", keytype = "GOALL",
                          multiVals = "list")
)

# Clean up results
go_symbol_list <- go_symbol_list[!is.na(go_symbol_list)]
go_symbol_list <- lapply(go_symbol_list, unique)

# Map gene symbols to Entrez IDs
entrez_id_map <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(unlist(go_symbol_list)),
                          column = "ENTREZID", keytype = "SYMBOL",
                          multiVals = "first")
)


##### Molecular Signatures Database (MSigDB) #########

# Need to ensure that the following collections are absent:
# C2 (curated), Canonical Pathways subcollection (Reactome subset retained)
# C4 (computational), Cancer Modules subcollection
# C5 (ontology), GO subcollection (redundant information)
# Don't want KEGG-containing parts of C2 or C4 for licensing reasons,
# and C5 contains GO (redundant)
colxn2keep  <- tolower(unique(msigdbr::msigdbr_collections()$gs_collection))
colxn2keep  <- colxn2keep[!colxn2keep %in% c("c5", "c9")] # Exclude entire C5 collection for now

# Some subcollections within the larger collections also need to be removed, for licensing reasons
subcolx2rmv <- c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:BIOCARTA", "CM") # CM = KEGG-derived cancer module

msigdb_df <- msigdbr::msigdbr(species = "Homo sapiens") |>
    mutate(group_code = tolower(gs_collection))

msigdb_filt <- msigdb_df |>
    dplyr::filter(group_code %in% colxn2keep) |>
    dplyr::filter(!gs_subcollection %in% subcolx2rmv) |>
    dplyr::select(gene_symbol, entrez_id = ncbi_gene,
                  pathway_id = gs_id, pathway_name = gs_name,
                  group_code, db_version)

msigdb_ver <- unique(msigdb_filt$db_version)

msigdb_filt <- msigdb_filt |>
    dplyr::select(-db_version) |> # No longer needed
    mutate(resource = "msigdb",
           collection = case_when(group_code == "h" ~ "hallmark",
                                  group_code == "c1" ~ "positional",
                                  group_code == "c2" ~ "curated",
                                  group_code == "c3" ~ "regulatory_target",
                                  group_code == "c4" ~ "computational",
                                  group_code == "c6" ~ "oncogenic_signature",
                                  group_code == "c7" ~ "immunologic_signature",
                                  group_code == "c8" ~ "cell_type_signature")) |>
    mutate(label = paste0("MSigDB ", toupper(group_code),
                          " (", stringr::str_to_title(gsub("\\_", " ", collection)), ")"))


# Confirm that the restricted subcollections aren't present
subcol_remaining <- table(msigdb_df[match(msigdb_filt$pathway_id, msigdb_df$gs_id), ]$gs_subcollection)
any(grepl("KEGG|CM", names(subcol_remaining)))


##############################
#--- Combine data frames ----#
##############################

# Must first convert GO data into a data frame
go_df <- stack(go_symbol_list)
go_df <- go_df |>
    dplyr::rename(gene_symbol = values,
                  pathway_id = ind) |>
    mutate(entrez_id = unname(entrez_id_map[gene_symbol]),
           pathway_name = go_name_map$TERM[match(pathway_id, go_name_map$GOID)],
           resource = "go",
           group_code = tolower(go_term_map$aspect[match(pathway_id, go_term_map$term)])) |>
    mutate(collection = dplyr::case_when(group_code == "bp" ~ "biological_process",
                                         group_code == "mf" ~ "molecular_function")) |>
    mutate(label = paste("GO", stringr::str_to_title(gsub("\\_", " ", collection))))

table(go_df$group_code)
head(go_df$label)
any(is.na(go_df$entrez_id))

# Making sure colnames match before combining
colnames(go_df)
colnames(msigdb_filt)

pathway_map <- dplyr::bind_rows(go_df, msigdb_filt) |>
    dplyr::select(gene_symbol, entrez_id, contains("pathway"), resource, collection, group_code, label)

head(pathway_map)
table(pathway_map$collection)
table(pathway_map$group_code)
table(pathway_map$label)


##############################
#--- Internal lookup data ---#
##############################

# Final pathway_map object is very large, with a lot of redundant information.
# This section will slim it down and remove non-unique mapping info that can
# be saved elsewhere.

group_code_lookup <- pathway_map |>
    dplyr::select(resource, collection, group_code, label) |>
    dplyr::distinct()

# One row per unique pathway: pathway_id  pathway_name
# Much smaller than embedding pathway_name in every gene row
path_name_lookup <- pathway_map |>
    dplyr::select(pathway_id, pathway_name) |>
    dplyr::distinct()

usethis::use_data(group_code_lookup, 
                  path_name_lookup, 
                  internal = TRUE, overwrite = TRUE)


##############################
#--- Reduce object size -----#
##############################

# Slim pathway_map to 4 columns; removed columns are reconstructable via
# above lookup objects.
# group_code is cast to a factor to reduce memory footprint further.
pathway_map <- pathway_map |>
  dplyr::select(gene_symbol, entrez_id, pathway_id, group_code) |>
  dplyr::mutate(group_code = factor(group_code))


##############################
#--- Save objects -----------#
##############################

# Don't forget to update download date & version info in R/data.R when
# re-saving this object!
save(pathway_map, file = "data/pathway_map.rda", compress = "gzip")

packageVersion("GO.db") # 3.23.0
packageVersion("msigdbr") # 26.1.1
msigdb_ver
