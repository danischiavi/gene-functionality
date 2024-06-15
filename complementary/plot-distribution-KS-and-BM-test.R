#### KOLMOROV-SMIRNOV TEST ANALYSIS ####

# PLOT TO VISUALIZE THE DISTRIBUTION OF FUNCTIONAL AND NEGATIVE-CONTROL FOR EACH FEATURE 
library(ggplot2)
library(dplyr)

# Load the data
# file_path <- "/Users/danischiavinato/Desktop/features/protein-features"
df <- read.csv(file_path, header=TRUE)

# Ensure 'Functional' is treated as a factor for plotting
df$Functional <- as.factor(df$Functional)

# List of feature names
feature_names <- setdiff(names(df), "Functional")

# Function to create plots for each feature
create_plots <- function(data, feature, target) {
  density_plot <- ggplot(data, aes(x = .data[[feature]], fill = .data[[target]])) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Density Plot of", feature, "by Functionality"),
         x = feature,
         y = "Density") +
    scale_fill_manual(values = c("blue", "red"), 
                      name = "Functionality", 
                      labels = c("No", "Yes")) +
    theme_minimal()
  
  box_plot <- ggplot(data, aes(x = .data[[target]], y = .data[[feature]], fill = .data[[target]])) +
    geom_boxplot() +
    labs(title = paste("Box Plot of", feature, "by Functionality"),
         x = "Functionality",
         y = feature) +
    scale_fill_manual(values = c("blue", "red"), 
                      name = "Functionality", 
                      labels = c("No", "Yes")) +
    theme_minimal()
  
  histogram_plot <- ggplot(data, aes(x = .data[[feature]], fill = .data[[target]])) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    facet_wrap(vars(.data[[target]]), scales = "free_y") +
    labs(title = paste("Histogram of", feature, "by Functionality"),
         x = feature,
         y = "Count") +
    scale_fill_manual(values = c("blue", "red"), 
                      name = "Functionality", 
                      labels = c("No", "Yes")) +
    theme_minimal()
  
  print(density_plot)
  print(box_plot)
  print(histogram_plot)
}

# Create plots for each feature
for (feature in feature_names) {
  create_plots(df, feature, "Functional")
}


#### SIGNED KOLMOGOROV-SMIRNOV TEST ####

# Load the data if not done yet
# file_path <- "/Users/danischiavinato/Desktop/features/protein-features"
df <- read.csv(file_path, header=TRUE)

# Define the 2 populations
functional_data <- data %>% filter(Functional == "Yes")
non_functional_data <- data %>% filter(Functional == "No")
# List of feature names (excluding 'Functional')
feature_names <- setdiff(names(data), "Functional")

# Function to calculate ECDF and difference
calculate_ecdf_difference <- function(feature) {
  functional_values <- functional_data[[feature]]
  non_functional_values <- non_functional_data[[feature]]
  
  ecdf_functional <- ecdf(functional_values)
  ecdf_non_functional <- ecdf(non_functional_values)
  
  data_range <- seq(min(c(functional_values, non_functional_values)),
                    max(c(functional_values, non_functional_values)),
                    length.out = 1000)
  
  cdf_functional <- ecdf_functional(data_range)
  cdf_non_functional <- ecdf_non_functional(data_range)
  
  cdf_difference <- cdf_functional - cdf_non_functional
  
  list(data_range = data_range, cdf_difference = cdf_difference)
}

# Function to perform permutation test
perform_permutation_test <- function(functional_values, non_functional_values, data_range, num_permutations = 1000) {
  observed_difference <- max(ecdf(functional_values)(data_range) - ecdf(non_functional_values)(data_range))
  
  permutation_differences <- numeric(num_permutations)
  
  combined_data <- c(functional_values, non_functional_values)
  n_functional <- length(functional_values)
  n_non_functional <- length(non_functional_values)
  
  for (i in 1:num_permutations) {
    shuffled_labels <- sample(rep(c("Functional", "Non-Functional"), times = c(n_functional, n_non_functional)))
    shuffled_functional_data <- combined_data[shuffled_labels == "Functional"]
    shuffled_non_functional_data <- combined_data[shuffled_labels == "Non-Functional"]
    
    shuffled_cdf_functional <- ecdf(shuffled_functional_data)
    shuffled_cdf_non_functional <- ecdf(shuffled_non_functional_data)
    
    shuffled_cdf_functional_values <- shuffled_cdf_functional(data_range)
    shuffled_cdf_non_functional_values <- shuffled_cdf_non_functional(data_range)
    
    shuffled_cdf_difference <- shuffled_cdf_functional_values - shuffled_cdf_non_functional_values
    
    permutation_differences[i] <- max(shuffled_cdf_difference)
  }
  
  p_value <- mean(permutation_differences >= observed_difference)
  p_value
}

# Iterate over each feature and calculate the ECDF differences and p-values
results <- list()
for (feature in feature_names) {
  ecdf_diff <- calculate_ecdf_difference(feature)
  p_value <- perform_permutation_test(functional_data[[feature]], non_functional_data[[feature]], ecdf_diff$data_range)
  
  results[[feature]] <- list(
    data_range = ecdf_diff$data_range,
    cdf_difference = ecdf_diff$cdf_difference,
    p_value = p_value
  )
  
  # Plot the difference between the CDFs
  plot(ecdf_diff$data_range, ecdf_diff$cdf_difference, type = "l", 
       xlab = feature, ylab = "Difference in CDFs",
       main = paste("Difference between CDFs for", feature))
  
  cat("Feature:", feature, "P-value:", p_value, "\n")
}

#### KS-test ###
# Initialize a results data frame
ks_test_results <- data.frame(Feature = character(), Statistic = numeric(), P.Value = character(), stringsAsFactors = FALSE)

# Loop through each feature
for (feature in features) {
  if (feature != "Functional") {
    # Remove rows with NA values for the current feature
    functional_feature_data <- functional_data[[feature]] %>% na.omit()
    non_functional_feature_data <- non_functional_data[[feature]] %>% na.omit()
   
    ks_test <- ks.test(functional_feature_data, non_functional_feature_data)
    # Append results to the data frame
    ks_test_results <- rbind(ks_test_results, data.frame(Feature = feature, Statistic = ks_test$statistic, P.Value = ks_test$p.value))
    }
    # Skip "Functional" feature
    next
}

print(ks_test_results)



data <- read.csv(file_path, header=TRUE)

# Ensure 'Functional' is treated as a factor for plotting
data$Functional <- as.factor(data$Functional)

# Define the 2 populations
functional_data <- data %>% filter(Functional == "Yes")
non_functional_data <- data %>% filter(Functional == "No")

# List of feature names (excluding 'Functional')
feature_names <- setdiff(names(data), "Functional")

# Function to calculate signed KS statistic
calculate_signed_ks_statistic <- function(feature) {
  functional_values <- functional_data[[feature]]
  non_functional_values <- non_functional_data[[feature]]
  
  ecdf_functional <- ecdf(functional_values)
  ecdf_non_functional <- ecdf(non_functional_values)
  
  data_range <- sort(unique(c(functional_values, non_functional_values)))
  
  cdf_functional <- ecdf_functional(data_range)
  cdf_non_functional <- ecdf_non_functional(data_range)
  
  cdf_difference <- cdf_functional - cdf_non_functional
  
  max_diff <- max(cdf_difference)
  min_diff <- min(cdf_difference)
  
  signed_D <- 0
  if(max_diff > abs(min_diff)) {
    signed_D <- max_diff
  } else {
    signed_D <- min_diff
  }
  
  list(data_range = data_range, cdf_difference = cdf_difference, signed_ks_statistic = signed_D)
}

# Function to perform permutation test
perform_permutation_test <- function(signed_D, functional_values, non_functional_values, data_range, num_permutations = 1000) {
  observed_difference <- signed_D
  
  permutation_differences <- numeric(num_permutations)
  
  combined_data <- c(functional_values, non_functional_values)
  n_functional <- length(functional_values)
  n_non_functional <- length(non_functional_values)
  
  for (i in 1:num_permutations) {
    shuffled_indices <- sample(length(combined_data))
    functional_values <- functional_data[["phyloP_max_214w"]]
    non_functional_values <- non_functional_data[["phyloP_max_214w"]]
    # Split the shuffled indices into functional and non-functional parts
    shuffled_functional_indices <- shuffled_indices[1:length(functional_values)]
    shuffled_non_functional_indices <- shuffled_indices[(length(functional_values) + 1):length(combined_data)]
    
    # Use shuffled indices to reorder the combined data
    shuffled_functional_data <- combined_data[shuffled_functional_indices]
    shuffled_non_functional_data <- combined_data[shuffled_non_functional_indices]
    shuffled_labels <- sample(rep(c("Functional", "Non-Functional"), times = c(n_functional, n_non_functional)))
    shuffled_functional_data <- combined_data[shuffled_labels == "Functional"]
    shuffled_non_functional_data <- combined_data[shuffled_labels == "Non-Functional"]
    
    shuffled_ecdf_functional <- ecdf(shuffled_functional_data)
    shuffled_ecdf_non_functional <- ecdf(shuffled_non_functional_data)
    
    shuffled_cdf_functional_values <- shuffled_ecdf_functional(data_range)
    shuffled_cdf_non_functional_values <- shuffled_ecdf_non_functional(data_range)
    
    shuffled_cdf_difference <- shuffled_cdf_functional_values - shuffled_cdf_non_functional_values
    
    shuffled_max_diff <- max(shuffled_cdf_difference)
    shuffled_min_diff <- min(shuffled_cdf_difference)
    
    shuffled_signed_D <- 0
    if(shuffled_max_diff > abs(shuffled_min_diff)) {
      shuffled_signed_D <- shuffled_max_diff
    } else {
      shuffled_signed_D <- shuffled_min_diff
    }
    
    permutation_differences[i] <- shuffled_signed_D
  }
  
  p_value <- mean(permutation_differences >= observed_difference)
  p_value
}

# Iterate over each feature and calculate the signed KS statistics and p-values
results <- list()
for (feature in feature_names) {
  ks_result <- calculate_signed_ks_statistic(feature)
  p_value <- perform_permutation_test(ks_result$signed_ks_statistic, functional_data[[feature]], non_functional_data[[feature]], ks_result$data_range)
  
  results[[feature]] <- list(
    data_range = ks_result$data_range,
    cdf_difference = ks_result$cdf_difference,
    signed_ks_statistic = ks_result$signed_ks_statistic,
    p_value = p_value
  )
  
  # Plot the difference between the CDFs
  #plot(ks_result$data_range, ks_result$cdf_difference, type = "l", 
   #    xlab = feature, ylab = "Difference in CDFs",
    #   main = paste("Difference between CDFs for", feature))
  
  cat("Feature:", feature, "Signed KS Statistic:", ks_result$signed_ks_statistic, "P-value:", p_value, "\n")
}
