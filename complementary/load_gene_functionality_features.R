load_gene_functionality_features <- function() {
  source("complementary/config.R")
  library(tidyverse)

  # Load existing features
  funcProtExon2Data <- data.frame(read.csv(FUNC_PROT_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID,-Strand)
  funcProtExon3Data <- data.frame(read.csv(FUNC_PROT_EXON3_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID,-Strand)
  funcLncrnaExon1Data <- data.frame(read.csv(FUNC_LNCRNA_EXON1_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID,-Strand)
  funcLncrnaExon2Data <- data.frame(read.csv(FUNC_LNCRNA_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID,-Strand)
  funcSncrnaData <- data.frame(read.csv(FUNC_SNCRNA_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-GeneID,-Strand)

  protExon2NCData <- data.frame(read.csv(NC_PROT_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene,-GeneID,-Strand)
  protExon3NCData <- data.frame(read.csv(NC_PROT_EXON3_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene,-GeneID,-Strand)
  lncrnaExon1NCData <- data.frame(read.csv(NC_LNCRNA_EXON1_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene,-GeneID,-Strand)
  lncrnaExon2NCData <- data.frame(read.csv(NC_LNCRNA_EXON2_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene,-GeneID,-Strand)
  sncrnaNCData <- data.frame(read.csv(NC_SNCRNA_FEATURES_FILE, header=TRUE)) %>% dplyr::select(-ID,-Functional,-Chromosome,-Start,-End,-Sequence,-DistanceGene,-GeneID,-Strand)

  feature_matrix <- rbind(data.frame(Dataset="protein-coding-exon2", funcProtExon2Data), 
                       data.frame(Dataset="protein-coding-exon3", funcProtExon3Data),
                       data.frame(Dataset = "lncrna-exon1", funcLncrnaExon1Data), 
                       data.frame(Dataset = "lncrna-exon2", funcLncrnaExon2Data),
                       data.frame(Dataset = "short-ncrna", funcSncrnaData),
                       
                       data.frame(Dataset = "protein-exon2-negative-control", protExon2NCData), 
                       data.frame(Dataset = "protein-exon3-negative-control", protExon3NCData),
                       data.frame(Dataset = "lncrna-exon1-negative-control", lncrnaExon1NCData), 
                       data.frame(Dataset = "lncrna-exon2-negative-control", lncrnaExon2NCData),
                       data.frame(Dataset = "short-ncrna-negative-control", sncrnaNCData)
  )

  return(feature_matrix)
}