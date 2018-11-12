library(GenomicRanges)
i <- 1
ov_info <- readRDS("../../data/ov_cell_line_info.rds")
sample <- ov_info[i,]
segDir <- "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/wgs_segments"
pgdx.id <- sample$PGDx.ID
seg.file <- file.path(segDir, paste0(pgdx.id, ".rds"))
segs <- readRDS(seg.file)

test <- sapply(1:(length(segs) - 1), function(i){
    m1 <- segs$seg.mean[i]
    m2 <- segs$seg.mean[i+ 1]
    return(m2 - m1)
})