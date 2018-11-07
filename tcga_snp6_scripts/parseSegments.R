library(dplyr)
# Read in all segments
segments <- read.table("../data/pancan12_absolute.segtab.txt", header = TRUE)
# Get segments for cases we care about
sample.info <- read.table("/dcl01/scharpf1/data/bams/TCGA/paths/tumor-normal-index",
                          header = TRUE)
samples <- as.character(sample.info$tumor.id)
all.samples <- substr(as.character(segments$Sample), 1, 16)
samples <- samples[samples %in% all.samples]
segments$Sample <- all.samples
for(i in 1:length(samples)){
    sample <- samples[i]
    sample.segs <- subset(segments, Sample == sample)
}

