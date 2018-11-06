library("jfkit")
project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
tcga.gbm.path <- file.path(project.path, "tcga_gbm")
normalizeToReference(reference.path = tcga.gbm.path, 
                     create.new.reference = TRUE,
                     count.flag.threshold = 0.2,
                     normalize.controls = TRUE)
quit('no')
