#install.packages("UpSetR")
#install.packages("dplyr") # For data manipulation, if needed
library(ggplot2)
library(dplyr)
library(UpSetR)

source("complementary/config.R")

# BED Files
db1_data <- read.csv("data/raw/gencodev44_chr22.bed", header=TRUE, sep = "\t")
db2_data <- read.csv("data/raw/gencodev45_chr22.bed", header=TRUE, sep = "\t")
db3_data <- read.csv("data/raw/hgnc_chr22.bed", header=TRUE, sep = "\t")
db4_data <- read.csv("data/raw/ncbi_chr22.bed", header=TRUE, sep = "\t")
db5_data <- read.csv("data/raw/rnacentral_chr22.bed", header=TRUE, sep = "\t")
db6_data <- read.csv("data/raw/uniprot_chr22.bed", header=TRUE, sep = "\t")


all_annotations <- rbind(db1_data, db2_data, db3_data, db4_data, db5_data, db6_data)
summary(all_annotations)
# Cast type from Integer to Numeric
all_annotations$Start <- as.numeric(all_annotations$Start)
all_annotations$End <- as.numeric(all_annotations$End)
class(all_annotations$Start)

# Rename NCBI to RefSeq
if ("NCBI" %in% names(all_annotations)) {
  names(all_annotations)[names(all_annotations) == 'NCBI'] <- 'RefSeq'
} else {
  warning("Column 'NCBI' not found — skipping rename. Actual columns: ", paste(names(all_annotations), collapse=", "))
}

cat("Column names after rename:\n")
print(names(all_annotations))

expected_sets <- c("GencodeV44","GencodeV45","HGNC","RefSeq","RNACentral","UniProt")
available_sets <- intersect(expected_sets, names(all_annotations))
missing_sets   <- setdiff(expected_sets, names(all_annotations))

if (length(missing_sets) > 0) {
  stop("The following set columns are missing from the data: ",
       paste(missing_sets, collapse=", "),
       "\nAvailable columns: ", paste(names(all_annotations), collapse=", "))
}

uplot <- upset(all_annotations, sets = available_sets, order.by = "freq",
               text.scale = c(3, 3, 2.5, 2, 3, 2.5))

uplot$Base <- uplot$Base + scale_x_continuous(n.breaks = 4, expand = expansion(mult = c(0, 0.25)))

# upset() uses grid graphics
# Equivalent to: scale=3, width=3840, height=2160, units="px", dpi=600
png(filename = file.path(VIOLIN_PLOTS_DIR, "png", "figureS1-Upset-plot.png"),
    width = 3840 * 3, height = 2160 * 3, res = 600, bg = "white")
print(uplot)
dev.off()

pdf(file = file.path(VIOLIN_PLOTS_DIR, "pdf", "figureS1-Upset-plot.pdf"),
    width = (3840 * 3) / 600, height = (2160 * 3) / 600)
print(uplot)
dev.off()
