library(jfkit)
library(rtracklayer)

sureselectv2 <- import("/dcl01/scharpf1/data/bams/TCGA/capture/SureSelect_All_Exon_V2_with_annotation.hg38.bed",
                     format = "bed")
sureselectv2 <- reduce(sureselectv2)
data(filters.hg38)
seqlevels(filters.hg38) <- paste0("chr", c(1:22, "X", "Y"))
seqlevels(sureselectv2) <- seqlevels(filters.hg38)
seqlengths(sureselectv2) <- seqlengths(filters.hg38)
genome(sureselectv2) <- "hg38"
sureselectv2 <- sureselectv2[width(sureselectv2) >= 25]

project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
if(!dir.exists(project.path)){
    dir.create(project.path, recursive = TRUE)
}
tcga.ov.path <- file.path(project.path, "tcga_ov")
if(!dir.exists(tcga.ov.path)){
    dir.create(tcga.ov.path, recursive = TRUE)
}

make_bins(targeted.ranges = sureselectv2,
          filters.hg38,
          reference.genome="hg38",
          insert.size = 600L, 
          reference.path = tcga.ov.path,
          min.gc = 0.3,
          max.gc = 0.7)
quit("no")
