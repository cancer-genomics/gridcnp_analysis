#### Note that I've only adopted the wholegenome plot to JFKit so far, not the individual gene plot
library(ggplot2)
library(stringr)
library(svfilters.hg18)

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
segDir <-  "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/jfkit/tumor_segments"
seg.files <- list.files(segDir)
### Choose sample
seg.file <- seg.files[i]

### Now same for bins
binDir <-  "/dcl01/scharpf1/data/gridcnp_analysis/ov_cell_lines/jfkit/cases"
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

cns <- readRDS(file.path(segDir, seg.file))
cns <- keepSeqlevels(cns, c(paste0("chr", seq(1,22,1)), "chrX", "chrY"), pruning.mode = "coarse")

#Making data.frame to plot   

cns <- data.frame(chromosome = as.character(seqnames(cns)), 
                  start = start(ranges(cns)),  
                  end = end(ranges(cns)),  
                  log2 = cns$seg.mean)
cns$chromosome <- factor(cns$chromosome, levels = c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))

purity <- "N/A"
sample <- str_extract(seg.file, "PGDX[0-9]{1,}T")


p <-
    ggplot() + 
    geom_point(data = cnr, 
               aes(x = (start+end)/2, y = log2, color = type), 
               alpha = 0.4, size = 0.3) +
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
    ggtitle(label = paste(sample, ": segmented log2R (Purity = ", purity, ") --JFKit/segmentr", sep = "")) +
    theme(plot.title = element_text(face = c("italic", "bold"), family = "Times", colour = "saddlebrown"))
plot.dir <- "../../visualizations/ov_cell_lines_segments"
if(!dir.exists(plot.dir)){
    dir.create(plot.dir, recursive = TRUE)
}
ggsave(file.path(plot.dir, paste0(sample, ".jpg")),
       plot = p, width = 11, height = 8.5, unit = "in")
quit('no')
