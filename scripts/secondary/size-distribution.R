# Calculates the size distribution of the exons 
lengths <- read.csv("/Users/danischiavinato/Desktop/gene-functionality/others/length")

z_scores <- scale(lengths$lncrna)
outliers <- abs(z_scores) > 3  # Adjust the threshold as needed
outlier_values <- lengths$lncrna[outliers]
cleaned_data <- lengths$lncrna[!outliers]
quantile(cleaned_data, probs = 0.9)
quantile(cleaned_data, probs = 0.1)
hist(cleaned_data, breaks = 1000, main = "exon (1&2) length distribution - lncrna", xlab = "length")
boxplot(cleaned_data,main = "exon (1&2) length distribution - lncrna") 

short <- na.omit(lengths$short)
z_scores <- scale(short)
outliers <- abs(z_scores) > 3  # Adjust the threshold as needed
outlier_values <- short[outliers]
cleaned_data <- short[!outliers]
quantile(cleaned_data, probs = 0.9)
quantile(cleaned_data, probs = 0.1)
hist(filtered_data, breaks = 1000, main = "length distribution - short-ncrna", xlab = "length")
boxplot(cleaned_data,main = "length distribution - short-ncrna") 
percentile_99 <- quantile(short, 0.99)

# Keep only values within the 99th percentile
filtered_data <- short[short <= percentile_99]

protein <- na.omit(lengths$protein)
percentile_99 <- quantile(protein, 0.99)
percentile_0.01 <- quantile(protein, 0.01)
filtered_data <- protein[protein <= percentile_99]
filtered_data <- filtered_data[filtered_data >= percentile_0.01]
quantile(filtered_data, probs = 0.9)
quantile(filtered_data, probs = 0.1)

hist(filtered_data_2, breaks = 1000, main = "exons (2&3) length distribution - protein-coding-rna", xlab = "length")
boxplot(filtered_data_2,main = "length distribution - protein-coding-rna") 