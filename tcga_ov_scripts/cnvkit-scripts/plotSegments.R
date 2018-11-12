library(ggplot2)
library(dplyr)

i <- as.numeric(commandArgs(trailingOnly = TRUE))
options(bitmapType='cairo')

#----------------------------------
# Getting locations of centromeres
#------------------------------------------------------------------------------------------------------------------
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
plot.info.centromeres <- data.frame(chromosome = seqnames(centromeres),
                                    start = start(ranges(centromeres)), 
                                    end = end(ranges(centromeres)))
plot.info.centromeres$chromosome <- as.character(plot.info.centromeres$chromosome)
plot.info.centromeres$chromosome <- factor(plot.info.centromeres$chromosome, levels = c(paste0("chr", seq(1,22,1)), "chrX", "chrY"))

### List files in segDir
segDir <- "/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/tumorsegs"
seg.files <- list.files(segDir)
### Choose sample
seg.file <- seg.files[i]

### Now same for bins
binDir <- "/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/tumornorm"
bin.file <- file.path(binDir, gsub(".cns", ".cnr", seg.file))

cnr <- read.table(bin.file, header = TRUE)

#Making data.frame to plot
cnr$chromosome <- factor(cnr$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))

### WXS segments
cns <- read.table(file.path(segDir, seg.file), header = TRUE)
#Making data.frame to plot   
cns$chromosome <- factor(cns$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))
cns$type <- rep("CNVkit", nrow(cns))
cns <- select(cns, -gene, -depth, -probes, -weight)

### WGS segments
wgsSegDir <- "/dcl01/scharpf1/data/gridcnp_analysis/tcga_wgs_segs"
cns.wgs <- readRDS(file.path(wgsSegDir, gsub(".cns", ".rds", seg.file)))
cns.wgs <- keepSeqlevels(cns.wgs, c(paste0("chr", seq(1,22,1)), "chrX", "chrY"), pruning.mode = "coarse")


cns.wgs <- data.frame(chromosome = as.character(seqnames(cns.wgs)), 
                      start = start(ranges(cns.wgs)),  
                      end = end(ranges(cns.wgs)),  
                      log2 = cns.wgs$seg.mean)
cns.wgs$chromosome <- factor(cns.wgs$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))
cns.wgs$type <- rep("SNP6", nrow(cns.wgs))
cns <- rbind(cns, cns.wgs)

sample <- gsub(".cns", "", seg.file)

p <-
    ggplot() + 
    geom_point(data = cnr, 
               aes(x = (start+end)/2, y = log2), 
               alpha = 0.2, size = 0.3) +
    geom_segment(data = cns, 
                 aes(x = start, xend = end, y = log2, yend = log2, color = type)) +
    geom_rect(data = plot.info.centromeres, 
              aes(xmin = start, xmax = end, ymin = -3, ymax = 3), 
              fill = "pink", alpha = 0.8) + 
    facet_wrap(~chromosome, scales = "free_x", nrow = 1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    ylab("log2R") + 
    xlab("Coordinate") + 
    ylim(-2.5, 2.5) +
    ggtitle(label = paste(sample, ": segmented log2R--CNVkit", sep = "")) +
    theme(plot.title = element_text(face = c("italic", "bold"), family = "Times", colour = "saddlebrown"))
plot.dir <- "../../visualizations/tcga_ov_segments_cnvkit"
if(!dir.exists(plot.dir)){
    dir.create(plot.dir, recursive = TRUE)
}
ggsave(file.path(plot.dir, paste0(sample, ".jpg")),
       plot = p, width = 11, height = 8.5, unit = "in")
quit('no')
