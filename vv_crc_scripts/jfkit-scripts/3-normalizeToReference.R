library("jfkit")
project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
vv.crc.path <- file.path(project.path, "vv_crc")
normalizeToReference(reference.path = vv.crc.path, 
                     create.new.reference = TRUE,
                     count.flag.threshold = 0.2,
                     normalize.controls = FALSE)
quit('no')
