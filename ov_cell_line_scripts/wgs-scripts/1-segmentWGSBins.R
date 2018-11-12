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

# Subsample to reduce hyper segmentation
set.seed(123)
chrs <- paste0("chr", c(1:22, "X", "Y"))
subsamps <- lapply(chrs, function(chr){
    bin.chr <- bins1kb[seqnames(bins1kb) == chr]
    subsamp <- sort(sample(1:length(bin.chr), round(1 / 10 * length(bin.chr)), replace = F))
    return(bin.chr[subsamp])
})
bins.subsamp <- unlist(GRangesList(subsamps))

segments <- segmentr::segmentBins(bins.subsamp, alpha = 0.001, undo.splits = "sdundo",
                                  undo.SD = 3, centromeres = centromeres)
seg.dir <- "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/wgs_segments"
if(!dir.exists(seg.dir)){
    dir.create(seg.dir, recursive = TRUE)
}
pgdx.id <- sample$PGDx.ID
seg.file <- file.path(seg.dir, paste0(pgdx.id, ".rds"))
saveRDS(segments, seg.file)
quit("no")