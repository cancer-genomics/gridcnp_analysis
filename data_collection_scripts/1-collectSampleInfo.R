##### This is a script to build .csv files with all info necessary
##### to 1) Download files 2) Perform analysis 
##### We want: POD location, cluster location, sample ID (lab, pgdx, tcga),
#####          probe information (type & file location), matched WGS info, project, matched normal info,
#####          batch information & whether we have a panel of normals done using same probe set
##### This script will read in various sources of sample info, and parse them
##### to build multiple .csv file with necessary pieces of info for each data type
##### We will have 3 .csv files: Ovarian cell lines, colorectal cancer (used in TEC-Seq paper), and TCGA

library(readxl)
library(tidyverse)

#### First Ovarian cell line data
#### Will read in data from Rob w/ lab IDs, then use master data to match that to PGDX IDs 
ov_cl_cp_seq <- read.csv("../data/amrf_091513_with_ids.csv", header = T)
#### Subset to unique values of lab id
ov_cl_cp_seq <- ov_cl_cp_seq[!duplicated(ov_cl_cp_seq$Lab.ID),]
ov_master_data <-
    read_excel("../data/Ovarian Sample Master Sheet 9.26.16.xlsx", skip = 1) %>%
    rename(Lab.ID = `Lab ID`, PGDx.ID = `PGDx ID`) %>%
    select(Lab.ID, PGDx.ID)
ov_cl_cp_seq <- inner_join(ov_cl_cp_seq, ov_master_data)
#### Now directory for certain PGDXids
#### Do blackslashes show up?
pgdx_613_658_targeted_dir <- "//10.113.118.248/data1/Data\\ from\\ PGDx/Ovarian/PGDX607-658\\ and\\ 1449"
ov_cell_line_file_data <-
    ov_cl_cp_seq %>%
    select(Lab.ID, PGDx.ID, Genome.Build) %>%
    filter(grepl("^PGDX6.*", PGDx.ID)) %>%
    mutate(targetd_seq_pod_location = file.path(pgdx_613_658_targeted_dir,
                                                paste0("t_", PGDx.ID, "_Cp.bam"))) %>%
    ### Should find out if this is CpPa or something different
    mutate(type = "ovarian", targeted_panel = "Cp", wgs_bin_loc = NA)
if(FALSE){
    ## Find path to bin-level data for above samples on cluster
    library(ovarian.manuscript)
    data(ovarian_ids)
    write_csv(ovarian_ids, "../data/ovarian_cell_line_manuscript.csv")
}
ovarian_ids <- read_csv("../data/ovarian_cell_line_manuscript.csv")
wgs.info <- ovarian_ids %>%
    right_join(ov_cell_line_file_data, by=c("alternative_id"="Lab.ID")) %>%
    mutate(bindir="/dcl01/scharpf/data/rscharpf/projects/OvarianCellLines/preprocess/3background_adj",
           wgs_bin_loc=file.path(bindir, paste0(alternative_id, ".bam.rds")))
stopifnot(identical(wgs.info$alternative_id, ov_cell_line_file_data$Lab.ID))
ov_cell_line_file_data <-
    ov_cell_line_file_data %>%
    mutate(wgs_bin_loc=wgs.info$wgs_bin_loc) %>%
    filter(file.exists(wgs_bin_loc))  ## drops one sample
saveRDS(ov_cell_line_file_data, file.path("..", "data", "ov_cell_line_info.rds"))
write.csv(ov_cell_line_file_data, file.path("..", "data", "ov_cell_line_info.csv"))
### Need panel of healthys that were done using same targeted sequencing
### Need WGS 

##############################################
### CRC samples
crc.sample.info <- read_excel("../data/091118_Cancer_Cases_with_Matched_T_N_P_WGS_and_Targeted_Seq.xlsx")
crc.targeted.dir <- "//10.113.118.248/data2/Phallen_et_al_EGA/"
crc.wgs.dir <- "/dcl01/scharpf1/data/bams/colorectal/tissue/wholegenome"
normal.wgs.dir <- "/dcl01/scharpf1/data/bams/healthy/tissue/wholegenome"
cluster.targeted.dir <- "/dcl01/scharpf1/data/bams/gridcnp_analysis/vv_crc"
### Remove first row (merged cells) and select only CRC cases
crc.sample.info <-
    crc.sample.info[-1,] %>%
    filter(`Patient Type` == "CRC") %>%
    select(Patient, `Patient Type`, Stage, starts_with("Copy"), starts_with("Fold"),
           X__1) %>%
    setNames(tolower(names(.))) %>%
    rename(labid = patient, cancer_type = `patient type`, cnv_genes = starts_with("copy"),
           cnv_fold = starts_with("fold"), pgdxid = x__1) %>%
    mutate(pgdxid = gsub("P$", "", pgdxid)) %>%
    mutate(targeted_tumor_pod_location = paste0(crc.targeted.dir, paste0("t_", pgdxid, "T_CpPa.bam"))) %>%
    mutate(targeted_normal_pod_location = gsub("/t_", "/n_", targeted_tumor_pod_location)) %>%
    mutate(wgs_tumor_cluster_location = file.path(crc.wgs.dir, paste0(pgdxid, "T_WGS_eland_sortednofa.bam"))) %>%
    mutate(wgs_normal_cluster_location = file.path(normal.wgs.dir, paste0(pgdxid, "N_WGS_eland_sortednofa.bam"))) %>%
    mutate(cluster_targeted_location = cluster.targeted.dir)
write.csv(crc.sample.info, "../data/crc_sample_info.csv")