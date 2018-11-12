library(dplyr)
library(readxl)
#### Read in capture kit data
tcga.ck <-
    read.csv("../data/TCGA_kit_metadata.csv", header = T) %>%
    rename(id = TCGA.ID, aliquot = Aliquot.ID, type = Tumor.Type, kit = Capture.Kit, date = Sequencing.Date) %>%
    mutate(date = substr(date, 1, 10)) %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"))
tcga.ck$id <- as.character(tcga.ck$id)
tcga.ck$type <- as.character(tcga.ck$type)
### subset by SureSelect Human All Exon 38 Mb v2 kit
tcga.ck <- filter(tcga.ck, kit == "SureSelect Human All Exon 38 Mb v2")

### Read in purity & ploidy
tcga.ploidy <- read.table("../data/pancan12.sample_info.txt", header = TRUE, fill = TRUE)
tcga.ploidy$tcga_id <- as.character(tcga.ploidy$tcga_id)
tcga.ploidy$id <- substr(tcga.ploidy$tcga_id, 1, 12)
tcga.ploidy <- rename(tcga.ploidy, type = disease)

## Merge, keeping maximum date, sequenced in 2011 or later
tcga.merged <-
    inner_join(tcga.ck, tcga.ploidy, by = c("id", "type")) %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(id) %>%
    filter(date == max(date)) %>%
    ungroup() %>%
    filter(year >= 2011) %>%
    arrange(type) 

### Subset to ovarian, as these are the most common type
ovarian <- filter(tcga.merged, type == "OV")

### Save date
write.csv(ovarian, "../data/potential_tcga_ovarian_samples.csv", row.names = FALSE)

