library(SomaDataIO)

# Example 11K (v5.0) adat from SomaLogic-Data GitHub repo
src_url <- "https://github.com/SomaLogic/SomaLogic-Data/raw/main/example_data_v5.0_plasma.adat"
tmp <- tempfile(fileext = ".adat")

download.file(url = src_url, destfile = tmp, quiet = TRUE)
example_data_11k <- SomaDataIO::read_adat(tmp)

# Save as global package object
save(example_data_11k, file = "data/example_data_11k.rda", compress = "xz")

unlink(tmp) # Cleanup
