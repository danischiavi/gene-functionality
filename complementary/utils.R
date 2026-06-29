# Utility functions for data manipulation and plotting

# Reads just the identifier column from each source file (in canonical Dataset
# order) and concatenates it, so it can be re-attached to a feature/z-score
# dataframe positionally without disturbing the loader functions that other
# scripts (heatmaps, PCA, KS stats) depend on returning purely numeric columns.
load_rnaids <- function(file_paths, id_column, sep = ",") {
  unlist(lapply(file_paths, function(path) {
    if (sep == "\t") read.delim(path, header = TRUE, sep = "\t")[[id_column]]
    else read.csv(path, header = TRUE)[[id_column]]
  }), use.names = FALSE)
}

melt_zscore_data <- function(dataframe, selected_features) {
  plot_list <- list()
  # Create an empty list to store data frames
  long_data_frames <- list()

  for (col in selected_features) {
    # Melt z-scores for violin plot
    df_with_zscores <- dataframe %>%
      pivot_longer(cols = col, names_to = "feature", values_to = "z_score") %>%
      dplyr::select(feature, z_score)

    # Store the result in the list
    long_data_frames[[col]] <- df_with_zscores
  }

  # Combine all melted data frames using bind_rows
  df_long <- dplyr::bind_rows(long_data_frames)
  return(df_long)
}

create_jitter_plot <- function(combined_df, plot_file_path) {
  p <- ggplot(combined_df, aes(x = feature, y = z_score, fill = Dataset)) +
    geom_jitter(aes(color = Dataset), width = 0.2, alpha = 0.5) +
    scale_color_manual(values = c(
      "mRNA (+)" = "blue", "lncRNA (+)" = "red", "sncRNA (+)" = "green",
      "mRNA (-)" = "lightblue", "lncRNA (-)" = "pink", "sncRNA (-)" = "lightgreen"
    )) +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(
      title = "Z-scores of Epigenetic Features",
      x = "Epigenetic Feature",
      y = "Robust Z-score"
    )

  ggsave(plot_file_path, p, width = 15, height = 10)
}

create_violin_plot <- function(df_long_format, ylimit, title, violin_adjust = 1, violin_scale = "width", plot_path) {
  p <- ggplot(df_long_format, aes(x = feature, y = z_score, fill = `Gene type`)) +
    #geom_jitter() + 
    geom_violin(scale = violin_scale, na.rm = TRUE, adjust = violin_adjust, trim = FALSE) +
    geom_boxplot(alpha=0.0, outliers=TRUE, na.rm = TRUE, position = position_dodge(width = 0.9), width=0.2) +
    facet_wrap(~ feature, scales = "free") + 
    labs(title = "Conservation", x = "Feature", y = "Robust Z-score") +
    theme_minimal(base_size = 34) +
    theme(
      axis.text.x = element_text(size = 0, angle = 45, hjust = 1),  # x-axis tick text hidden; facet strip is the label
      axis.text.y = element_text(size = 30),  # y-axis tick text size
      axis.title.x = element_text(size = 0),  # x-axis title hidden (redundant with facet strip)
      axis.title.y = element_text(size = 38),  # y-axis title size
      legend.position = "none",
      legend.title = element_text(size = 20),  # legend title size
      legend.text = element_text(size = 18),  # legend text size
      plot.title = element_text(size = 0, hjust = 0.5), # plot title hidden (facet strip is the label)
      plot.subtitle = element_blank(),
      axis.line.x = element_blank(),          # Remove x-axis line
      axis.ticks.x = element_blank(),         # Remove x-axis ticks
      panel.grid.major.x = element_blank(),   # Remove major grid lines along x-axis
      panel.grid.minor.x = element_blank(),    # Remove minor grid lines along x-axis
      strip.text = element_text(size = 40, margin = ggplot2::margin(0,0,40,0))
    ) +
    scale_fill_manual(values = c("#F4A582FF","#c9e3f6FF","#D6604DFF","#56bdfcFF","#e37b88FF","#53a4f5FF")) + # colors
    coord_cartesian(ylim = ylimit) # Zoom into an specific area without removing data points
  
  ggsave(plot_path, p, scale = 3, width = 3840, height = 2160, units = "px", bg = "white", dpi = 600)
  return(p)
}

################################
# Function to perform K-S test #
run_ks_tests <- function(dataN, dataP) {
  n <- ncol(dataN)
  results <- matrix(nrow = 4, ncol = n)
  colnames(results) <- colnames(dataN)
  rownames(results) <- c("signed_D", "max", "min", "p.val")
  
  for (i in colnames(dataN)) {
    positive_col <- dataP[[i]]
    negative_col <- dataN[[i]]
    
    ks_test <- custom_ks_test(negative_col, positive_col)
    ks_test_p <- ks.test(negative_col, positive_col)
    
    results[1, i] <- ks_test$signed_D
    results[2, i] <- ks_test$max_diff
    results[3, i] <- ks_test$min_diff
    results[4, i] <- ks_test_p$p.value
  }
  return(results)
}

############################
# Custom K-S test function #
custom_ks_test <- function(x, y) {
  # Sort data
  x <- sort(x)
  y <- sort(y)
  
  # Calculate empirical cumulative distribution functions (ECDF)
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  
  # Get the unique values from both samples
  unique_vals <- sort(unique(c(x, y)))
  
  # Calculate the raw differences between ECDFs
  diffs <- ecdf_x(unique_vals) - ecdf_y(unique_vals)
  
  # Find the maximum difference (positive or negative)
  max_diff <- max(diffs)
  min_diff <- min(diffs)
  
  signed_D <- if (max_diff > abs(min_diff)) max_diff else min_diff
  
  # Return the maximum and minimum differences
  return(list(signed_D = signed_D, max_diff = max_diff, min_diff = min_diff))
}
