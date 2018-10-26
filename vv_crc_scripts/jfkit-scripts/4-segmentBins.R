i <- as.numeric(commandArgs(trailingOnly = TRUE))
library(svfilters.hg19)
library(segmentr)

data(gaps, package = "svfilters.hg19")
centromeres = gaps[which(gaps$type == "centromere")]
inDir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/jfkit/cases"
### Currently 66 of these files
files <- list.files(inDir, "*log2norm.rds")
file <- files[i]
bins <- readRDS(file.path(inDir, file))

# Removing flagged bins 
bins <- bins[bins$flag == FALSE]

# Remove bins with variance above 95th percentile in normals
upper.var <- quantile(bins$var, .95)
bins <- bins[bins$var < upper.var]

set.seed(123)
segments <- segmentr::segmentBins(bins,alpha = 0.01, undo.splits = "sdundo",
                                  undo.SD = 2, centromeres = centromeres)

# Not fine-tuning segments with n.probes = 3
hits <- which(segments$n.probes <= 3)
if (length(hits) > 0) {
   lte3 <- segments[hits]
   segments <- segments[-hits]
} else {
   lte3 <- NULL
}

finetuned <- segmentr::finetune.segments(bins = bins, segments = segments,
                                         alpha = 0.01, centromeres = centromeres, undo.SD = 1)

if (length(lte3) > 0) {
  finetuned <- sort(c(finetuned, lte3))
}

output.file <- gsub("log2norm", "", file)
segDir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/jfkit/tumor_segments"
if(!dir.exists(segDir)) {
    dir.create(segDir, recursive = TRUE)
}
saveRDS(finetuned, file = file.path(segDir, output.file))
quit('no')
