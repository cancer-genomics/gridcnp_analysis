library(trellis)
library(GenomicAlignments)
bamfile <- commandArgs(trailingOnly = TRUE)

data(bins1kb, package="svfilters.hg19")
bins <- bins1kb
flags <- scanBamFlag(isDuplicate=FALSE,
                     isPaired=TRUE,
                     isUnmappedQuery=FALSE,
                     hasUnmappedMate=FALSE,
                     isProperPair=TRUE)
bviews <- BamViews(bamPaths=bamfile, bamRanges=bins)
bins$cnt <- binnedCounts(bviews)
bins <- bins[ bins$cnt > 0 ]
bins$std_cnt <- binNormalize(bins)
set.seed(123)
bins$log_ratio <- binGCCorrect(bins)
datadir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/wgs/bincounts"
if(!dir.exists(datadir)){
    dir.create(datadir, recursive = TRUE)
}
binfile <- file.path(datadir, gsub(".bam", ".rds", basename(bamfile)))
saveRDS(bins, binfile)
quit("no")
