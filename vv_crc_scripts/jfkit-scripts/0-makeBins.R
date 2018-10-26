library(jfkit)

CpPa2_hg_19_ranges <- readRDS("/dcl01/scharpf1/data/probe_sets/probe_granges/CpPa2_reduced_hg19.rds")
data(filters.hg19)

project.path <- "/dcl01/scharpf1/data/gridcnp_analysis"
if(!dir.exists(project.path)){
    dir.create(project.path, recursive = TRUE)
}
vv.crc.path <- file.path(project.path, "vv_crc")
if(!dir.exists(vv.crc.path)){
    dir.create(vv.crc.path, recursive = TRUE)
}

make_bins(targeted.ranges = CpPa2_hg_19_ranges,
          filters.hg19,
          reference.genome="hg19",
          insert.size = 600L, 
          reference.path = vv.crc.path,
          min.gc = 0.3,
          max.gc = 0.7)
quit("no")
