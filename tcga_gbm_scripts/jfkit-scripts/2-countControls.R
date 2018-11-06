library(jfkit)
library(getopt)
library(readxl)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
id <- options.args[1]
bamdir <- options.args[2]
bam.file <- file.path(bamdir, id)

############
### See if we've already completed this run for JFkit
project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
tcga.gbm.path <- file.path(project.path, "tcga_gbm")
jfkit.dir <- file.path(tcga.gbm.path, "jfkit/controls")
jfkit.file <- file.path(jfkit.dir, gsub(".bam", "log2.rds", id))
jfkit.exists <- file.exists(jfkit.file)
if(jfkit.exists) {
    quit("no")
}
#############
countExperiment(bam.file, 
                is.case = FALSE, 
                reference.path = tcga.gbm.path)
quit("no")
