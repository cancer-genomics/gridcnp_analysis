library(GenomicRanges)

control1 <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls/TCGA-06-0184-10Blog2.rds")
control2 <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls/TCGA-12-1089-10Clog2.rds")


sum(control1$rawcounts[control1$type == "targeted"]) / sum(control1$rawcounts)
sum(control2$rawcounts[control2$type == "targeted"]) / sum(control2$rawcounts)

sum(control1$rawcounts[control1$type == "targeted"] == 0)
sum(control2$rawcounts[control2$type == "targeted"] == 0)

plot(control1$adjusted[control1$type == "targeted"]- control2$adjusted[control2$type == "targeted"])

plot(control1$adjusted[control1$type == "off.target"] - control2$adjusted[control2$type == "off.target"])

targeted.bins <- control1$type=="targeted"
control1 <- control1[targeted.bins]
control2 <- control2[targeted.bins]

plot(density(control1$weightedcounts))
lines(density(control2$weightedcounts), col = "red")

plot(density(width(control1)))
plot(density((width(control1[control1$rawcounts==0]))))

plot(density(width(control2)))
plot(density((width(control2[control2$rawcounts==0]))))


autosomal <- seqnames(control1) %in% paste0("chr", 1:22)
sum(control1[autosomal]$rawcounts == 0)
sum(control2[autosomal]$rawcounts == 0)

b <- control1[control1$rawcounts > 500 & control2$rawcounts == 0]
par(mfrow=c(2,1))
plot(control1$adjusted)
plot(control2$adjusted)

### Let's look at the density of raw counts for all normals

normals <- list.files("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls/",
                      "log2.rds$", full.names = TRUE)

bin.list <- lapply(normals, function(f) {
    b <- readRDS(f)
    t <- b[b$type == "targeted"]
    df <- data.frame(bin = 1:length(t), adjcount = t$rawcounts / sum(t$rawcounts))
    df$sample <- rep(gsub(".rds", "", basename(f)), nrow(df))
    return(df)
})

library(ggplot2)
bin.df <- do.call(rbind, bin.list)
messy.samps <- c("TCGA-12-1089", "TCGA-12-3651", "TCGA-19-1788", "TCGA-19-2629")
q <- quantile(bin.df$adjcount, .95)

unique.samps <- unique(bin.df$sample)
is.messy <- sapply(unique.samps, function(s) {
    sub <- substr(s, 1, 12)
    return(sum(grepl(sub, messy.samps) > 0))
})
m <- match(bin.df$sample, names(is.messy))

bin.df$is.messy <- unname(is.messy)[m]
bin.df$sample2 <- ifelse(bin.df$is.messy == 1, bin.df$sample, "Not messy and/or no case")

### four samples with clear outliers
tcga.gbm.norms <-
    ggplot(bin.df, aes(adjcount, color = sample2, group = sample)) +
    geom_density() +
    xlim(0, q)
ggsave("../../visualizations/tcga_gbm_normal_count_density.jpg",
       plot = tcga.gbm.norms,
       width = 12, height = 8)

### What about if we do this with GBM vs OV

gbm.normals <- list.files("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls/",
                        "log2.rds$", full.names = TRUE) 
ov.normals <- list.files("/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/jfkit/controls/",
                          "log2.rds$", full.names = TRUE) 

all.normals <- c(gbm.normals, ov.normals)
bin.list <- lapply(all.normals, function(f) {
    b <- readRDS(f)
    t <- b[b$type == "targeted"]
    df <- data.frame(bin = 1:length(t), adjcount = t$rawcounts / sum(t$rawcounts))
    df$sample <- rep(gsub(".rds", "", basename(f)), nrow(df))
    df$type <- ifelse(f %in% gbm.normals, "GBM", "OV")
    return(df)
})
bin.df <- do.call(rbind, bin.list)
bin.df <-
    bin.df %>%
    group_by(sample) %>%
    mutate(nzero = sum(adjcount == 0))
bin.df <- subset(bin.df, nzero < 50000)
x.max <- quantile(bin.df$adjcount, .95)
gbm.vs.ov <-
    ggplot(bin.df, aes(adjcount, color = type, group = sample)) +
    geom_density() +
    xlim(0, q)
ggsave("../../visualizations/tcga_ov_vs_gbm_normal_count_density.jpg",
       plot = gbm.vs.ov,
       width = 12, height = 8)

n1  <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls//TCGA-12-3651-10Alog2.rds")
n2 <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls//TCGA-19-2629-10Alog2.rds")
n3 <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls//TCGA-12-1089-10Clog2.rds")
n4 <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/tcga_gbm/jfkit/controls//TCGA-06-0122-10Alog2.rds")

plot(n1$adjusted - n2$adjusted)
plot(n1$adjusted - n3$adjusted)
plot(n1$adjusted - n4$adjusted)

sum(n1$rawcounts[n1$type=="targeted"])/sum(n1$rawcounts)
sum(n2$rawcounts[n2$type=="targeted"])/sum(n2$rawcounts)
sum(n3$rawcounts[n3$type=="targeted"])/sum(n3$rawcounts)
sum(n4$rawcounts[n4$type=="targeted"])/sum(n4$rawcounts)


### Is this a problem with bed file?
suppressMessages(library(GenomicRanges))
file <- '/dcl01/scharpf1/data/bams/TCGA/capture/SureSelect_All_Exon_V2_with_annotation.hg19.minimal.bed'
a = read.table(file, header = FALSE, sep = '\t')
g <- GRanges(seqnames = a$V1, ranges = IRanges(start = a$V2, end = a$V3))


### What does this look like for targeted sequencing
crc.bins <- readRDS("/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/jfkit/controls/n_PGDX5881T_CpPa_processedlog2.rds")
min(crc.bins$rawcounts[crc.bins$type=="targeted"])
plot(density(crc.bins$weightedcounts[crc.bins$type=="targeted"]))
