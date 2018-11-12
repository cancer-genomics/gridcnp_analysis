#### Note that I've only adopted the wholegenome plot to JFKit so far, not the individual gene plot
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
centromeres <- keepSeqlevels(centromeres, c(paste0("chr", seq(1,22,1)), "chrX"), pruning.mode = "coarse")
plot.info.centromeres <- data.frame(chromosome = seqnames(centromeres),
                                    start = start(ranges(centromeres)), 
                                    end = end(ranges(centromeres)))
plot.info.centromeres$chromosome <- as.character(plot.info.centromeres$chromosome)
plot.info.centromeres$chromosome <- factor(plot.info.centromeres$chromosome, levels = c(paste0("chr", seq(1,22,1)), "chrX"))

### Get bins and segments for sample
### Read in WGS bin counts
ov_info <- readRDS("../../data/ov_cell_line_info.rds")
sample <- ov_info[i,]
wgs.bins <- readRDS(sample$wgs_bin_loc)
### Get actual bins
data(bins1kb)
bins1kb <- keepSeqlevels(bins1kb, paste0("chr", c(1:22, "X")), pruning.mode = "coarse")
### Have to divide by 1000
bins1kb$adjusted <- wgs.bins / 1000


segDir <- "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/wgs_segments"
pgdx.id <- sample$PGDx.ID
seg.file <- file.path(segDir, paste0(pgdx.id, ".rds"))

### Make cnr
cnr <- bins1kb

#Making data.frame to plot
cnr <- data.frame(chromosome = as.character(seqnames(cnr)),
                  start = start(ranges(cnr)), 
                  end = end(ranges(cnr)), 
                  log2 = cnr$adjusted)
cnr$chromosome <- factor(cnr$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX"))

cns <- readRDS(seg.file)
cns <- keepSeqlevels(cns, c(paste0("chr", seq(1,22,1)), "chrX"), pruning.mode = "coarse")

#Making data.frame to plot   

cns <- data.frame(chromosome = as.character(seqnames(cns)), 
                  start = start(ranges(cns)),  
                  end = end(ranges(cns)),  
                  log2 = cns$seg.mean)

cns$chromosome <- factor(cns$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX"))

#purity <- "N/A"
sample <- str_extract(seg.file, "PGDX[0-9]{1,}T")
#cnr <- cnr[seq(1, nrow(cnr), by = 10),]

p <-
    ggplot() + 
    geom_point(data = cnr, 
               aes(x = (start+end)/2, y = log2), 
               alpha = 0.1, size = 0.2, color = "grey66") +
    geom_segment(data = cns, 
                 aes(x = start, xend = end, y = log2, yend = log2), 
                 color = "black") +
    geom_rect(data = plot.info.centromeres, 
              aes(xmin = start, xmax = end, ymin = -3, ymax = 3), 
              fill = "pink", alpha = 0.8) + 
    facet_wrap(~chromosome, scales = "free_x", nrow = 1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    ylab("log2R") + 
    xlab("Coordinate") + 
    ylim(-3, 3) +
    ggtitle(label = sample) +
    theme(plot.title = element_text(face = c("italic", "bold"), family = "Times", colour = "saddlebrown"))
plot.dir <- "../../visualizations/ov_cell_lines_wgs_segments"
if(!dir.exists(plot.dir)){
    dir.create(plot.dir, recursive = TRUE)
}
ggsave(file.path(plot.dir, paste0(sample, ".jpg")),
       plot = p, width = 11, height = 8.5, unit = "in")
quit('no')
