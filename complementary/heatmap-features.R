#### FEATURES CORRELATION MATRIX - FEATURE:FEATURE (for each type of RNA individually) ####

# file <- "/Users/danischiavinato/Desktop/data-from-server/20240320/short-features"

data <- read.csv(file, header = TRUE)

data$Functional <- factor(data$Functional, levels=c("No","Yes")) 
functional_numeric <- as.numeric(factor(data$Functional)) - 1 
data$Functional <- functional_numeric
data <- na.roughfix(data)

cor_matrix <- cor(data, method = "spearman")

# Heatmap
library(ggplot2)
library(reshape2)

# Melt the correlation matrix into long format
cor_melted <- melt(cor_matrix)

ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Spearman Correlation") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    plot.title = element_text(size = 12),  # Increase title size
    axis.title = element_text(size = 10),  # Increase axis title size
    axis.text = element_text(size = 20),    # Increase axis label size
    legend.title = element_text(size = 10),  # Increase legend title size
    legend.text = element_text(size = 10),   # Increase legend text size
    legend.key.size = unit(2, "lines")       # Increase legend key size
  ) +
  labs(title = "Features Spearman Correlation - short-ncRNA", x = "Features", y = "Features")