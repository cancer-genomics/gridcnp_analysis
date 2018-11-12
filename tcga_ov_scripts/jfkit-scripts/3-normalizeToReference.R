library("jfkit")
project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
tcga.ov.path <- file.path(project.path, "tcga_ov")
normalizeToReference(reference.path = tcga.ov.path, 
                     create.new.reference = TRUE,
                     count.flag.threshold = 0.2,
                     normalize.controls = TRUE)
quit('no')
