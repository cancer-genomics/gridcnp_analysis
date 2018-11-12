library(dplyr)
library(GenomicRanges)
library(rtracklayer)
# Read in all segments
segments <- read.table("/dcl01/scharpf1/data/gridcnp_analysis/pancan12_relative.seg.txt", header = TRUE)
segments$Chromosome <- as.character(segments$Chromosome)
test <-
    segments %>%
    group_by(Sample)%>%
    summarise(nchrom = length(unique(Chromosome)),
              x = sum(Chromosome == "X"))
segments$Chromosome[segments$Chromosome == "23"] <- "Y"
segments$Sample <- substr(as.character(segments$Sample), 1, 16)
# Get segments for cases we care about
sample.info <- read.table("/dcl01/scharpf1/data/bams/TCGA/paths/tumor-normal-index",
                          header = TRUE)
samples <- as.character(sample.info$tumor.id)
samples <- samples[samples %in% unique(segments$Sample)]
seg.dir <- "/dcl01/scharpf1/data/gridcnp_analysis/tcga_wgs_segs"
if(!dir.exists(seg.dir)){
    dir.create(seg.dir, recursive = TRUE)
}
hg19_to_hg38 <- import.chain("/dcl01/scharpf1/data/liftover_chains/hg19ToHg38.over.chain")


for(i in 1:length(samples)){
    sample <- samples[i]
    sample.segs <- subset(segments, Sample == sample)
    segs.granges <- GRanges(seqnames = paste0("chr", sample.segs$Chromosome),
                            ranges = IRanges(start = sample.segs$Start, end = sample.segs$End),
                            nprobes = sample.segs$Num_Probes, seg.mean = sample.segs$Segment_Mean)
    segs.granges$id <- 1:length(segs.granges)
    ### Liftover to hg38
    segs.granges.hg38 <- liftOver(segs.granges, hg19_to_hg38)
    ### Keep most common chromosome for each segment
    segs.granges.hg38 <- lapply(segs.granges.hg38, function(ranges) {
        t <- table(seqnames(ranges))
        seq.max <- names(t[which.max(t)])
        ranges <- ranges[seqnames(ranges) == seq.max]
        return(ranges)
    })
    segs.granges.hg38 <- unlist(GRangesList(segs.granges.hg38))
    segs.granges.hg38.df <- data.frame(chr = seqnames(segs.granges.hg38),
                                       start = start(segs.granges.hg38),
                                       end = end(segs.granges.hg38),
                                       id = segs.granges.hg38$id,
                                       seg.mean = segs.granges.hg38$seg.mean)
    segs.granges.hg38.df <-
        segs.granges.hg38.df %>%
        group_by(id) %>%
        summarise(chr = unique(chr),
                  start = min(start),
                  end = max(end),
                  seg.mean = unique(seg.mean))
    segs.granges.hg38 <- GRanges(seqnames = segs.granges.hg38.df$chr,
                                 ranges = IRanges(start = segs.granges.hg38.df$start,
                                                  end = segs.granges.hg38.df$end),
                                 seg.mean = segs.granges.hg38.df$seg.mean)
    seqlevels(segs.granges) <- paste0("chr", c(1:22, "X", "Y"))
    genome(segs.granges) <-"hg38"
    output.file <- file.path(seg.dir, paste0(sample, ".rds"))
    saveRDS(segs.granges, output.file)
}
quit("no")

