library(dplyr)
#### Read in capture kit data
data <- read.csv("../data/TCGA_kit_metadata.csv", header = T)
data$TCGA.ID <- as.character(data$TCGA.ID)
data$Capture.Kit <- as.character(data$Capture.Kit)

### Get samples of interest
samples <- read.table("/dcl01/scharpf1/data/bams/TCGA/paths/tumor-normal-index",
                      header = T)
samples$tumor.id <- substr(samples$tumor.id, 1, 12)
colnames(samples)[1] <- "TCGA.ID"

data <- data[data$TCGA.ID %in% samples$TCGA.ID,]

merged.data <- inner_join(data, samples)
merged.data$Sequencing.Date <- as.character(merged.data$Sequencing.Date)
merged.data$Sequencing.Date <- substr(merged.data$Sequencing.Date, 1, 10)


messy.samps <- c("TCGA-12-1089", "TCGA-12-3651", "TCGA-19-1788", "TCGA-19-2629")
### Really bad ones are TCGA-12-1089 and TCGA-19-1788
gbm.data <-
    merged.data %>%
    filter(tumor.type == "GBM") %>%
    mutate(is.messy = TCGA.ID %in% messy.samps) %>%
    arrange(desc(Sequencing.Date))
gbm.data$dup <- duplicated(gbm.data$TCGA.ID)

write.csv(gbm.data, "../data/tcga_gbm_sequencing_data.csv", row.names = FALSE)

ov.data <- 
    merged.data %>%
    filter(tumor.type == "OV")
write.csv(ov.data, "../data/tcga_ov_sequencing_data.csv", row.names = FALSE)

