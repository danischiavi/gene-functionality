file <- "lncrna_epigenetic_features_matrix.csv"
file <- "short_ncrna_epigenetic_features_matrix.csv"
file <- "protein_epigenetic_features_matrix.csv"

data <- read.csv(file, header = TRUE, sep = "\t")
data$Functionality <- factor(data$Functionality, levels=c("No","Yes"))
functional_numeric <- as.numeric(factor(data$Functionality)) - 1
data$Functionality <- functional_numeric
data <- na.roughfix(data)
cor_matrix <- cor(data, method = "spearman")
# Heatmap
library(ggplot2)
library(reshape2)
# Melt the correlation matrix into long format
cor_melted <- melt(cor_matrix)


###################
ggplot(cor_melted, aes(Var1, Var2)) +  
  geom_tile(aes(fill = 1), color = "grey", show.legend = FALSE) + 
  geom_text(aes(Var1, Var2, label = round(value, 2), color = value), nudge_x = 0, nudge_y = 0, fontface = "bold") + 
  scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red", midpoint = 0, limits = c(-0.8, 0.8)) +
  scale_fill_gradient2(low = "grey", high = "grey", midpoint = 0) +  # For background
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    plot.title = element_text(size = 16),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 10),  
    legend.key.size = unit(2, "lines") 
  ) +
  labs(title = "Features Spearman Correlation - protein coding", x = "Features", y = "Features")
##################
