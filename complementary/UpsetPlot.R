install.packages("UpSetR")
install.packages("dplyr") # For data manipulation, if needed

# BED Files
db1_data <- read.csv("result/gencodev44.bed", header=TRUE, sep = "\t")
db2_data <- read.csv("result/gencodev45.bed", header=TRUE, sep = "\t")
db3_data <- read.csv("result/hgnc.bed", header=TRUE, sep = "\t")
db4_data <- read.csv("result/ncbi.bed", header=TRUE, sep = "\t")
db5_data <- read.csv("result/rnacentral.bed", header=TRUE, sep = "\t")
db6_data <- read.csv("result/uniprot.bed", header=TRUE, sep = "\t")

library(dplyr)
all_annotations <- rbind(db1_data, db2_data, db3_data, db4_data, db5_data, db6_data)
summary(all_annotations)
# Cast type from Integer to Numeric
all_annotations$Start <- as.numeric(all_annotations$Start)
all_annotations$End <- as.numeric(all_annotations$End)
class(all_annotations$Start)

library(UpSetR)
#help(upset)
#Join the first three columns into a new column named 'Region'
#all_annotations <- all_annotations %>% mutate(Region = paste(Chromosome, Start, End, sep=":"))

upset(all_annotations, sets = c("GencodeV44","GencodeV45","HGNC","NCBI","RNACentral","UniProt"), 
      order.by = "freq")
