i <- as.numeric(commandArgs(trailingOnly = TRUE))
library(segmentr)
library(rCGH)
library(GenomicRanges)

### get hg38 centromeres
hg38.centromeres <- rCGH::hg38
chr <- paste0("chr", hg38.centromeres$chrom) 
chr[chr=="chr23"] <- "chrX"
chr[chr=="chr24"] <- "chrY"
centromeres <- GRanges(seqnames = chr,
                       ranges = IRanges(start = hg38.centromeres$centromerStart - 100,
                                        end = hg38.centromeres$centromerEnd + 100))
inDir <- "/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/cases"
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
segDir <- "/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/tumor_segments"
if(!dir.exists(segDir)) {
    dir.create(segDir, recursive = TRUE)
}
saveRDS(finetuned, file = file.path(segDir, output.file))
quit('no')
