library(jfkit)

Cp_hg_18_ranges <- readRDS("/dcl01/scharpf1/data/probe_sets/probe_granges/Cp_reduced_hg18.rds")
data(filters.hg18)

project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
if(!dir.exists(project.path)){
    dir.create(project.path, recursive = TRUE)
}
ov.path <- file.path(project.path, "ov_cell_lines")
if(!dir.exists(ov.path)){
    dir.create(ov.path, recursive = TRUE)
}

make_bins(targeted.ranges = Cp_hg_18_ranges,
          filters.hg18,
          reference.genome="hg18",
          insert.size = 600L, 
          reference.path = ov.path,
          min.gc = 0.3,
          max.gc = 0.7)
quit("no")
