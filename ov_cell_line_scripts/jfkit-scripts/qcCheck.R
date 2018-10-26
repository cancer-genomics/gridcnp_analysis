library(GenomicRanges)
bin.dir <- "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/jfkit/cases"
files <- list.files(bin.dir, full.names = TRUE)
samples <- seq(from = 1L, to = 20L, by = 4)
sample.list <- lapply(samples, function(i) {
    file <- files[i]
    bins <- readRDS(file)
    df <- data.frame()
    
})
sample <- readRDS(list.files(bin.dir, full.names = TRUE)[1])