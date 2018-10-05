##### This is a script to build .csv files with all info necessary
##### to 1) Download files 2) Perform analysis 
##### We want: POD location, cluster location, sample ID (lab, pgdx, tcga),
#####          probe information (type & file location), matched WGS info, project, matched normal info,
#####          batch information & whether we have a panel of normals done using same probe set
##### This script will read in various sources of sample info, and parse them
##### to build multiple .csv file with necessary pieces of info for each data type
##### We will have 3 .csv files: Ovarian cell lines, colorectal cancer (used in TEC-Seq paper), and TCGA

library(readxl)
library(dplyr)

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
pgdx_613_658_targeted_dir <- "//10.113.118.248/data1/Data\ from\ PGDx/Ovarian/PGDX607-658\ and\ 1449"
ov_cell_line_file_data <-
    ov_cl_cp_seq %>%
    select(Lab.ID, PGDx.ID, Genome.Build) %>%
    filter(grepl("^PGDX6.*", PGDx.ID)) %>%
    mutate(targetd_seq_pod_location = file.path(pgdx_613_658_targeted_dir,
                                                paste0("t_", PGDx.ID, "_Cp.bam"))) %>%
    ### Should find out if this is CpPa or something different
    mutate(type = "ovarian", targeted_panel = "Cp")
### Need panel of healthys that were done using same targeted sequencing
### Need WGS 
    
