load_gene_functionality_zscores <- function() {
  source("complementary/config.R")
  library(dplyr)
  library(tidyr)

  # Load zscores files
  funcProtExon2Zscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "functional-protein-exon2-dataset-zscores.csv"), header=TRUE))
  funcProtExon3Zscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "functional-protein-exon3-dataset-zscores.csv"), header=TRUE))
  funcLncrnaExon1Zscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "functional-lncrna-exon1-dataset-zscores.csv"), header=TRUE))
  funcLncrnaExon2Zscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "functional-lncrna-exon2-dataset-zscores.csv"), header=TRUE))
  funcSncrnaZscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "functional-short-ncrna-dataset-zscores.csv"), header=TRUE))

  protExon2NCZscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "protein-exon2-negative-control-dataset-zscores.csv"), header=TRUE))
  protExon3NCZscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "protein-exon3-negative-control-dataset-zscores.csv"), header=TRUE))
  lncrnaExon1NCZscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "lncrna-exon1-negative-control-dataset-zscores.csv"), header=TRUE))
  lncrnaExon2NCZscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "lncrna-exon2-negative-control-dataset-zscores.csv"), header=TRUE))
  sncrnaNCZscores <- data.frame(read.csv(file.path(ZSCORES_DIR, "short-ncrna-negative-control-dataset-zscores.csv"), header=TRUE))

  # Join in a single dataframe
  zscores_all <- rbind(data.frame(Dataset="protein-coding-exon2", funcProtExon2Zscores), 
                        data.frame(Dataset="protein-coding-exon3", funcProtExon3Zscores),
                        data.frame(Dataset = "lncrna-exon1", funcLncrnaExon1Zscores), 
                        data.frame(Dataset = "lncrna-exon2", funcLncrnaExon2Zscores),
                        data.frame(Dataset = "short-ncrna", funcSncrnaZscores),
                        
                        data.frame(Dataset = "protein-exon2-negative-control", protExon2NCZscores), 
                        data.frame(Dataset = "protein-exon3-negative-control", protExon3NCZscores),
                        data.frame(Dataset = "lncrna-exon1-negative-control", lncrnaExon1NCZscores), 
                        data.frame(Dataset = "lncrna-exon2-negative-control", lncrnaExon2NCZscores),
                        data.frame(Dataset = "short-ncrna-negative-control", sncrnaNCZscores)
  )

  return(zscores_all)
}