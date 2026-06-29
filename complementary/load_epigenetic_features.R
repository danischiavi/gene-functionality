load_epigenetic_features <- function() {
  source("complementary/config.R")
  library(dplyr)
  library(tidyr)

  # Load existing histone features
  funcProtExon2HistoneData <- data.frame(read.csv(COMBINED_PROT_EXON2_HISTONE_FILE, header=TRUE))
  funcProtExon3HistoneData <- data.frame(read.csv(COMBINED_PROT_EXON3_HISTONE_FILE, header=TRUE))
  funcLncrnaExon1HistoneData <- data.frame(read.csv(COMBINED_LNCRNA_EXON1_HISTONE_FILE, header=TRUE))
  funcLncrnaExon2HistoneData <- data.frame(read.csv(COMBINED_LNCRNA_EXON2_HISTONE_FILE, header=TRUE))
  funcSncrnaHistoneData <- data.frame(read.csv(COMBINED_SNCRNA_HISTONE_FILE, header=TRUE))

  protExon2NCHistoneData <- data.frame(read.csv(COMBINED_NC_PROT_EXON2_HISTONE_FILE, header=TRUE))
  protExon3NCHistoneData <- data.frame(read.csv(COMBINED_NC_PROT_EXON3_HISTONE_FILE, header=TRUE))
  lncrnaExon1NCHistoneData <- data.frame(read.csv(COMBINED_NC_LNCRNA_EXON1_HISTONE_FILE, header=TRUE))
  lncrnaExon2NCHistoneData <- data.frame(read.csv(COMBINED_NC_LNCRNA_EXON2_HISTONE_FILE, header=TRUE))
  sncrnaNCHistoneData <- data.frame(read.csv(COMBINED_NC_SNCRNA_HISTONE_FILE, header=TRUE))

  # Load chromatin accessibility feature
  funcProtExon2chrmacc <- read.csv(CHR_ACC_PROT_EXON2_FILE, header = TRUE)
  funcProtExon3chrmacc <- read.csv(CHR_ACC_PROT_EXON3_FILE, header = TRUE)
  funcLncrnaExon1chrmacc <- read.csv(CHR_ACC_LNCRNA_EXON1_FILE, header = TRUE)
  funcLncrnaExon2chrmacc <- read.csv(CHR_ACC_LNCRNA_EXON2_FILE, header = TRUE)
  funcSncrnachrmacc <- read.csv(CHR_ACC_SNCRNA_FILE, header = TRUE)

  negProtExon2chrmacc <- read.csv(NC_CHR_ACC_PROT_EXON2_FILE, header = TRUE)
  negProtExon3chrmacc <- read.csv(NC_CHR_ACC_PROT_EXON3_FILE, header = TRUE)
  negLncrnaExon1chrmacc <- read.csv(NC_CHR_ACC_LNCRNA_EXON1_FILE, header = TRUE)
  negLncrnaExon2chrmacc <- read.csv(NC_CHR_ACC_LNCRNA_EXON2_FILE, header = TRUE)
  negSncrnachrmacc <- read.csv(NC_CHR_ACC_SNCRNA_FILE, header = TRUE)

  # Paste chromatin into epigenetic dataframes
  funcProtExon2HistoneData$chrm_acc_MaxScaledSignal <- funcProtExon2chrmacc$chrm_acc_MaxScaledSignal
  funcProtExon3HistoneData$chrm_acc_MaxScaledSignal <- funcProtExon3chrmacc$chrm_acc_MaxScaledSignal
  funcLncrnaExon1HistoneData$chrm_acc_MaxScaledSignal <- funcLncrnaExon1chrmacc$chrm_acc_MaxScaledSignal
  funcLncrnaExon2HistoneData$chrm_acc_MaxScaledSignal <- funcLncrnaExon2chrmacc$chrm_acc_MaxScaledSignal
  funcSncrnaHistoneData$chrm_acc_MaxScaledSignal <- funcSncrnachrmacc$chrm_acc_MaxScaledSignal

  protExon2NCHistoneData$chrm_acc_MaxScaledSignal <- negProtExon2chrmacc$chrm_acc_MaxScaledSignal
  protExon3NCHistoneData$chrm_acc_MaxScaledSignal <- negProtExon3chrmacc$chrm_acc_MaxScaledSignal
  lncrnaExon1NCHistoneData$chrm_acc_MaxScaledSignal <- negLncrnaExon1chrmacc$chrm_acc_MaxScaledSignal
  lncrnaExon2NCHistoneData$chrm_acc_MaxScaledSignal <- negLncrnaExon2chrmacc$chrm_acc_MaxScaledSignal
  sncrnaNCHistoneData$chrm_acc_MaxScaledSignal <- negSncrnachrmacc$chrm_acc_MaxScaledSignal

  # Load methylome feature
  funcProtExon2methylome <- read.csv(METHYLOME_PROT_EXON2_FILE, header = TRUE)
  funcProtExon3methylome <- read.csv(METHYLOME_PROT_EXON3_FILE, header = TRUE)
  funcLncrnaExon1methylome <- read.csv(METHYLOME_LNCRNA_EXON1_FILE, header = TRUE)
  funcLncrnaExon2methylome <- read.csv(METHYLOME_LNCRNA_EXON2_FILE, header = TRUE)
  funcSncrnamethylome <- read.csv(METHYLOME_SNCRNA_FILE, header = TRUE)

  negProtExon2methylome <- read.csv(NC_METHYLOME_PROT_EXON2_FILE, header = TRUE)
  negProtExon3methylome <- read.csv(NC_METHYLOME_PROT_EXON3_FILE, header = TRUE)
  negLncrnaExon1methylome <- read.csv(NC_METHYLOME_LNCRNA_EXON1_FILE, header = TRUE)
  negLncrnaExon2methylome <- read.csv(NC_METHYLOME_LNCRNA_EXON2_FILE, header = TRUE)
  negSncrnamethylome <- read.csv(NC_METHYLOME_SNCRNA_FILE, header = TRUE)

  # Integrate methylome into epigenetic dataframe
  funcProtExon2HistoneData$methylome <- funcProtExon2methylome$methylome
  funcProtExon3HistoneData$methylome <- funcProtExon3methylome$methylome
  funcLncrnaExon1HistoneData$methylome <- funcLncrnaExon1methylome$methylome
  funcLncrnaExon2HistoneData$methylome <- funcLncrnaExon2methylome$methylome
  funcSncrnaHistoneData$methylome <- funcSncrnamethylome$methylome

  protExon2NCHistoneData$methylome <- negProtExon2methylome$methylome
  protExon3NCHistoneData$methylome <- negProtExon3methylome$methylome
  lncrnaExon1NCHistoneData$methylome <- negLncrnaExon1methylome$methylome
  lncrnaExon2NCHistoneData$methylome <- negLncrnaExon2methylome$methylome
  sncrnaNCHistoneData$methylome <- negSncrnamethylome$methylome

  # Join all data in single csv file for convenience.
  cols1 <- colnames(funcProtExon2HistoneData)
  cols2 <- colnames(protExon2NCHistoneData)
  colsdiff1 <- setdiff(cols1, cols2)
  colsdiff2 <- setdiff(cols2, cols1)

  feature_matrix_epigenetic <- rbind(data.frame(Dataset="protein-coding-exon2", funcProtExon2HistoneData[, !(names(funcProtExon2HistoneData) %in% colsdiff1)]), 
                                    data.frame(Dataset="protein-coding-exon3", funcProtExon3HistoneData[, !(names(funcProtExon3HistoneData) %in% colsdiff1)]),
                                    data.frame(Dataset = "lncrna-exon1", funcLncrnaExon1HistoneData[, !(names(funcLncrnaExon1HistoneData) %in% colsdiff1)]), 
                                    data.frame(Dataset = "lncrna-exon2", funcLncrnaExon2HistoneData[, !(names(funcLncrnaExon2HistoneData) %in% colsdiff1)]),
                                    data.frame(Dataset = "short-ncrna", funcSncrnaHistoneData[, !(names(funcSncrnaHistoneData) %in% colsdiff1)]),
                                    
                                    data.frame(Dataset = "protein-exon2-negative-control", protExon2NCHistoneData[, !(names(protExon2NCHistoneData) %in% colsdiff2)]), 
                                    data.frame(Dataset = "protein-exon3-negative-control", protExon3NCHistoneData[, !(names(protExon3NCHistoneData) %in% colsdiff2)]),
                                    data.frame(Dataset = "lncrna-exon1-negative-control", lncrnaExon1NCHistoneData[, !(names(lncrnaExon1NCHistoneData) %in% colsdiff2)]), 
                                    data.frame(Dataset = "lncrna-exon2-negative-control", lncrnaExon2NCHistoneData[, !(names(lncrnaExon2NCHistoneData) %in% colsdiff2)]),
                                    data.frame(Dataset = "short-ncrna-negative-control", sncrnaNCHistoneData[, !(names(sncrnaNCHistoneData) %in% colsdiff2)])
  )

  return(feature_matrix_epigenetic)
}