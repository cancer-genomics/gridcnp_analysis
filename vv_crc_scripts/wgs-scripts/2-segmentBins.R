library(trellis)
library(GenomicRanges)
i <- as.numeric(commandArgs(trailingOnly = TRUE))
### Read in WGS bin counts
bin.files <- list.files("/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/wgs/bincounts",
                        full.names = TRUE) 
sample.file <- bin.files[i]
bins1kb <- readRDS(sample.file)
segs <- segmentBins(bins1kb, param = SegmentParam())
seg.dir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/wgs/segments"
if(!dir.exists(seg.dir)){
    dir.create(seg.dir, recursive = TRUE)
}
seg.file <- file.path(seg.dir, basename(sample.file))
saveRDS(segs, seg.file)
quit("no")