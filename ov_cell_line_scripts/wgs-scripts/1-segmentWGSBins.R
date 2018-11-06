library(segmentr)
library(svfilters.hg19)
library(GenomicRanges)
i <- as.numeric(commandArgs(trailingOnly = TRUE))
### Read in WGS bin counts
ov_info <- readRDS("../../data/ov_cell_line_info.rds")
sample <- ov_info[i,]
wgs.bins <- readRDS(sample$wgs_bin_loc)
### Get actual bins
data(bins1kb)
bins1kb <- keepSeqlevels(bins1kb, paste0("chr", c(1:22, "X")), pruning.mode = "coarse")
### Have to divide by 1000
bins1kb$adjusted <- wgs.bins / 1000

### Get centromeres
data(gaps, package = "svfilters.hg19")
centromeres = gaps[which(gaps$type == "centromere")]

### Drop chrY from centromeres
centromeres <- keepSeqlevels(centromeres, paste0("chr", c(1:22, "X")), pruning.mode = "coarse")

#bins1kb <- bins1kb[seq(1, length(bins1kb), by = 100)]
segments <- segmentr::segmentBins(bins1kb, alpha = 0.001, undo.splits = "sdundo",
                                  undo.SD = 5, centromeres = centromeres)
seg.dir <- "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/wgs_segments"
if(!dir.exists(seg.dir)){
    dir.create(seg.dir, recursive = TRUE)
}
pgdx.id <- sample$PGDx.ID
seg.file <- file.path(seg.dir, paste0(pgdx.id, ".rds"))
saveRDS(segments, seg.file)
quit("no")