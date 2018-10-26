library("jfkit")
project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
ov.path <- file.path(project.path, "ov_cell_lines")

normalizeToReference(reference.path = ov.path, 
                     create.new.reference = TRUE,
                     flat.reference = TRUE,
                     count.flag.threshold = 0.2,
                     normalize.controls = FALSE)
quit('no')
