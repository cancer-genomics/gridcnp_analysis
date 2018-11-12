library(ggplot2)
library(stringr)
library(svfilters.hg19)

i <- as.numeric(commandArgs(trailingOnly = TRUE))


options(bitmapType='cairo')

data(gaps)


#----------------------------------
# Getting locations of centromeres
#------------------------------------------------------------------------------------------------------------------
centromeres <- subset(gaps, gaps$type == "centromere")
centromeres <- keepSeqlevels(centromeres, c(paste0("chr", seq(1,22,1)), "chrX", "chrY"), pruning.mode = "coarse")
plot.info.centromeres <- data.frame(chromosome = seqnames(centromeres),
                                    start = start(ranges(centromeres)), 
                                    end = end(ranges(centromeres)))
plot.info.centromeres$chromosome <- as.character(plot.info.centromeres$chromosome)
plot.info.centromeres$chromosome <- factor(plot.info.centromeres$chromosome, levels = c(paste0("chr", seq(1,22,1)), "chrX", "chrY"))

### List files in segDir
segDir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/jfkit/tumor_segments"
seg.files <- list.files(segDir)
### Choose sample
seg.file <- seg.files[i]

### Now same for bins
binDir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/jfkit/cases"
bin.file <- file.path(binDir, gsub(".rds", "log2norm.rds", seg.file))

cnr <- readRDS(bin.file)
## Filter CNR like we did for segmenting
cnr <- cnr[cnr$flag == FALSE & cnr$var < 0.5]

#Making data.frame to plot

cnr <- data.frame(chromosome = as.character(seqnames(cnr)),
                  start = start(ranges(cnr)), 
                  end = end(ranges(cnr)), 
                  log2 = cnr$adjusted, 
                  type = cnr$type)
cnr$chromosome <- factor(cnr$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))

### Targeted segments

cns <- readRDS(file.path(segDir, seg.file))
cns <- keepSeqlevels(cns, c(paste0("chr", seq(1,22,1)), "chrX", "chrY"), pruning.mode = "coarse")

#Making data.frame to plot   

cns <- data.frame(chromosome = as.character(seqnames(cns)), 
                  start = start(ranges(cns)),  
                  end = end(ranges(cns)),  
                  log2 = cns$seg.mean)
cns$type <- rep("JFKit", nrow(cns))
### WGS segments
### WGS segments
wgsSegDir <- "/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/wgs/segments"
sample <- str_extract(seg.file, "PGDX[0-9]{1,}T")
cns.wgs <- readRDS(file.path(wgsSegDir, paste0(sample, "_WGS_eland_processed.rds")))
cns.wgs <- keepSeqlevels(cns.wgs, c(paste0("chr", seq(1,22,1)), "chrX", "chrY"), pruning.mode = "coarse")


cns.wgs <- data.frame(chromosome = as.character(seqnames(cns.wgs)), 
                      start = start(ranges(cns.wgs)),  
                      end = end(ranges(cns.wgs)),  
                      log2 = cns.wgs$seg.mean)
cns.wgs$chromosome <- factor(cns.wgs$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))
cns.wgs$type <- rep("WGS", nrow(cns.wgs))
cns <- rbind(cns, cns.wgs)
cns$chromosome <- factor(cns$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))

p <-
    ggplot() + 
    geom_point(data = cnr, 
               aes(x = (start+end)/2, y = log2), 
               alpha = 0.3, size = 0.2,
               color = "grey66") +
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
    ggtitle(label = paste(sample, ": segmented log2R --JFKit/segmentr", sep = "")) +
    theme(plot.title = element_text(face = c("italic", "bold"), family = "Times", colour = "saddlebrown"))
plot.dir <- "../../visualizations/vv_crc_segments"
if(!dir.exists(plot.dir)){
    dir.create(plot.dir, recursive = TRUE)
}
ggsave(file.path(plot.dir, paste0(sample, ".jpg")),
       plot = p, width = 11, height = 8.5, unit = "in")
quit('no')
