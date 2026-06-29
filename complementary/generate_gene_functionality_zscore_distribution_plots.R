# Create violin plots from gene functionality z-score data.
# Load necessary libraries.
library(ggplot2)
library(tidyr)
library(patchwork)
library(scales)

source("complementary/config.R")
source("complementary/load_gene_functionality_zscores.R")
source("complementary/utils.R")
zscores_all <- load_gene_functionality_zscores()

PLOT_SELECTED_FEATURES <- c("Random","GC_percentage",
                            "CpG","GA","GT","TA",
                            "phyloP_max_241w","phyloP_max_100w",
                            "GERP_91_mammals_max","GERP_63_amniotes_max",
                            "fickett","coding_potential","Max_covariance",
                            "MAF_avg",
                            "MFE","Interaction_ave",
                            "SNP_density",
                            "methylome",
                            "T","C","G","MaxEnt_splicing")

PLOT_SELECTED_FEATURES_LABELS <- c("Random number","GC%",
                                   "CpG","GA","GT","TA",
                                   "PhyloP-mammals","PhyloP-vertebrates",
                                   "GERP-mammals","GERP-vertebrates",
                                   "Fickett","RNAcode","Covariance",
                                   "MAF",
                                   "MFE","Interactions",
                                   "SNPs",
                                   "Methylome",
                                   "T%","C%","G%","MaxEnt splicing")


PLOT_SELECTED_LOG <- c("lowComplexity_density",
                       "RPKM_tissue","RPKM_primary.cell",
                       "copy_number","repeat_distance",
                       "H3K9ac_MaxScaledSignal","H3K79me1_MaxScaledSignal","H3K79me2_MaxScaledSignal","chrm_acc_MaxScaledSignal")

PLOT_SELECTED_LOG_LABELS <- c("Low Complexity",
                       "Tissue RPKM","Primary Cell RPKM",
                       "Copies","Repeat Free",
                       "H3K9ac","H3K79me1","H3K79me2","Chromatin Accessibility")

# Function to create distribution plots
create_distribution_plots <- function(data, gene_type, color, output_dir = "plots", axis_label = "Robust Z-score", type = "Functional", title, tag) {
  # Create output directories if they don't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "png"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "pdf"), showWarnings = FALSE, recursive = TRUE)

  # Prepare data for plotting
  long_data <- data %>%
    tidyr::pivot_longer(everything(), names_to = "Feature", values_to = "Value") %>%
    mutate(Type = type)

  # List of features to plot
  features <- unique(long_data$Feature)
  plots_list <- list()

  # Create plots for each feature
  for (feature in features) {
    feature_data <- long_data %>% filter(Feature == feature)
    # Remove NAs
    feature_data <- feature_data[!is.na(feature_data$Value), ]
    # Ensure column is numeric
    feature_data$Value <- as.numeric(feature_data$Value)

    # 1. Histogram
    p1 <- ggplot(feature_data, aes(x = Value, fill = Type, colour = Type)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
      labs(x = axis_label, y = "Count") +
      scale_fill_manual(values = c(color)) +
      scale_color_manual(values = c(color)) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
      #                   breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
      #                   labels = scales::label_number()
      #) +
      theme_minimal(base_size = 84) +
      theme(axis.text.x = element_text(size = 36, angle = 45),
            axis.text.y = element_text(size = 36),  # Increase y-axis text size
            axis.title.x = element_text(size = 44),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"
      )

    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_histogram.png")), p1, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)

    # 2. Density plot
    p2 <- ggplot(feature_data, aes(x = Value, fill = Type, colour = Type)) +
      geom_density(alpha = 0.6) +
      theme_minimal(base_size = 84) +
      labs(x = axis_label, y = "Density") +
      scale_fill_manual(values = c(color)) +
      scale_color_manual(values = c(color)) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      #scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
      #                   breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
      #                   labels = scales::label_number()
      #) +
      theme(axis.text.x = element_text(size = 36, angle = 45),
            axis.text.y = element_text(size = 36),  # Increase y-axis text size
            axis.title.x = element_text(size = 44),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"#,
            #plot.title = element_text(size = 64, face = "italic", hjust = 0.5)  # center the title
      )

    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_density.png")), p2, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)

    # 3. Boxplot
    p3 <- ggplot(feature_data, aes(x = feature, y = Value, fill = Type)) +
      geom_boxplot(outlier.colour = color, outlier.size = 5) +
      theme_minimal(base_size = 84) +
      labs(y = axis_label) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      #scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      #                   breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
      #                   labels = scales::label_number()
      #) +
      theme(axis.text.x = element_text(size = 0, angle = 45),
            axis.text.y = element_text(size = 36),  # Increase y-axis text size
            axis.title.x = element_text(size = 0),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"
      ) +
      scale_fill_manual(values = c(color))

    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_boxplot.png")), p3, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)

    # 4. Join all plots in one
    combined_plot <- (p1 + p2 + p3)

    combined_plot_annotated <- wrap_elements(combined_plot +
                                               plot_annotation(title = title,
                                                               tag_levels = list(tag)) &
                                               theme(plot.title = element_text(size = 64,
                                                                               face = "italic",
                                                                               hjust = 0.5),
                                                     strip.text = element_text(margin = margin(0,0,30,0)),
                                                     plot.tag.position = c(0, 1),
                                                     plot.tag = element_text(size = 38, face = "bold", hjust = 0, vjust = 0, margin = margin(40,40,40,40))))

    #ggsave(paste0(gene_type, "_", feature, "_combined.png"), combined_plot_annotated,
    #       path = file.path(output_dir, "png"), scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    #ggsave(paste0(gene_type, "_", feature, "_combined.pdf"), combined_plot_annotated,
    #       path = file.path(output_dir, "pdf"), scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    plots_list[[feature]] <- combined_plot_annotated
  }
  return(plots_list)
}

# Function to create log scale distribution plots
create_distribution_log_plots <- function(data, gene_type, color, output_dir = "plots", axis_label = "Robust Z-score", type = "Functional", title, tag) {
  # Create output directories if they don't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "png"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "pdf"), showWarnings = FALSE, recursive = TRUE)
  
  # Prepare data for plotting
  long_data <- data %>%
    tidyr::pivot_longer(everything(), names_to = "Feature", values_to = "Value") %>%
    mutate(Type = type)
  
  # List of features to plot
  features <- unique(long_data$Feature)
  plots_list <- list()
  
  # Create plots for each feature
  for (feature in features) {
    feature_data <- long_data %>% filter(Feature == feature)
    # Remove NAs
    feature_data <- feature_data[!is.na(feature_data$Value), ]
    # Ensure column is numeric
    feature_data$Value <- as.numeric(feature_data$Value)
    
    # 1. Histogram
    p1 <- ggplot(feature_data, aes(x = Value, fill = Type, colour = Type)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
      labs(x = axis_label, y = "Count") +
      scale_fill_manual(values = c(color)) +
      scale_color_manual(values = c(color)) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
                         labels = scales::label_number()
      ) +
      theme_minimal(base_size = 84) +
      theme(axis.text.x = element_text(size = 36, angle = 45),
            axis.text.y = element_text(size = 36),  # Increase y-axis text size
            axis.title.x = element_text(size = 44),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"
      )
    
    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_histogram.png")), p1, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    
    # 2. Density plot
    p2 <- ggplot(feature_data, aes(x = Value, fill = Type, colour = Type)) +
      geom_density(alpha = 0.6) +
      theme_minimal(base_size = 84) +
      labs(x = axis_label, y = "Density") +
      scale_fill_manual(values = c(color)) +
      scale_color_manual(values = c(color)) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
                         labels = scales::label_number()
      ) +
      theme(axis.text.x = element_text(size = 36, angle = 45),
            axis.text.y = element_text(size = 36),  # Increase y-axis text size
            axis.title.x = element_text(size = 44),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"#,
            #plot.title = element_text(size = 64, face = "italic", hjust = 0.5)  # center the title
      )
    
    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_density.png")), p2, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    
    # 3. Boxplot
    p3 <- ggplot(feature_data, aes(x = feature, y = Value, fill = Type)) +
      geom_boxplot(outlier.colour = color, outlier.size = 5) +
      theme_minimal(base_size = 84) +
      labs(y = axis_label) +
      #scale_x_log10(breaks = scales::log_breaks(n = 10)) +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(-1000, -100, -10, -1, 10, 100, 1000, 10000, 100000,1000000,10000000),
                         labels = scales::label_number()
      ) +
      theme(axis.text.x = element_text(size = 0, angle = 45),
            axis.text.y = element_text(size = 36),  # Increase y-axis text size
            axis.title.x = element_text(size = 0),  # Increase x-axis title size
            axis.title.y = element_text(size = 44),  # Increase y-axis title size
            legend.position = "none"
      ) +
      scale_fill_manual(values = c(color))
    
    #ggsave(file.path(output_dir, paste0(gene_type, "_", feature, "_boxplot.png")), p3, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    
    # 4. Join all plots in one
    combined_plot <- (p1 + p2 + p3)
    
    combined_plot_annotated <- wrap_elements(combined_plot + 
                                               plot_annotation(title = title,
                                                               tag_levels = list(tag)) &
                                               theme(plot.title = element_text(size = 64,
                                                                               face = "italic",
                                                                               hjust = 0.5),
                                                     strip.text = element_text(margin = margin(0,0,30,0)),
                                                     plot.tag.position = c(0, 1),
                                                     plot.tag = element_text(size = 38, face = "bold", hjust = 0, vjust = 0, margin = margin(40,40,40,40))))
    
    #ggsave(paste0(gene_type, "_", feature, "_combined.png"), combined_plot_annotated,
    #       path = file.path(output_dir, "png"), scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    #ggsave(paste0(gene_type, "_", feature, "_combined.pdf"), combined_plot_annotated,
    #       path = file.path(output_dir, "pdf"), scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
    plots_list[[feature]] <- combined_plot_annotated
  }
  return(plots_list)
}


# Generate plots to visualize distributions:
unique(zscores_all$Dataset)
p1 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_mrna", color = "#F4A582FF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR, title = "mRNA(+)", tag = c('A','B','C'))
p2 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_mrna_negative", "#c9e3f6FF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR,
                                type = "Negative Control", title = "mRNA(-)", tag = c('D','E','F'))

p3 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("short-ncrna")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_sncRNA", color = "#D6603FFF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR, title = "sncRNA(+)", tag = c('A','B','C'))
p4 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("short-ncrna-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_sncRNA_negative", color = "#56bdfcFF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR,
                                type = "Negative Control", title = "sncRNA(-)", tag = c('D','E','F'))

p5 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_lncRNA", color = "#e37b88FF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR, title = "lncRNA(+)", tag = c('A','B','C'))
p6 <- create_distribution_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_FEATURES)),
                                "zscores_lncRNA_negative", color = "#5271ffFF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR,
                                type = "Negative Control", title = "lncRNA(-)", tag = c('D','E','F'))


# Log scale plots
p1_log <- create_distribution_log_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("protein-coding-exon2", "protein-coding-exon3")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_LOG)),
                                "zscores_mrna", color = "#F4A582FF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR, title = "mRNA(+)", tag = c('A','B','C'))
p2_log <- create_distribution_log_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("protein-exon2-negative-control", "protein-exon3-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_LOG)),
                                "zscores_mrna_negative", "#c9e3f6FF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR,
                                type = "Negative Control", title = "mRNA(-)", tag = c('D','E','F'))

p3_log <- create_distribution_log_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("short-ncrna")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_LOG)),
                                "zscores_sncRNA", color = "#D6603FFF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR, title = "sncRNA(+)", tag = c('A','B','C'))
p4_log <- create_distribution_log_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("short-ncrna-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_LOG)),
                                "zscores_sncRNA_negative", color = "#56bdfcFF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR,
                                type = "Negative Control", title = "sncRNA(-)", tag = c('D','E','F'))

p5_log <- create_distribution_log_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("lncrna-exon1", "lncrna-exon2")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_LOG)),
                                "zscores_lncRNA", color = "#e37b88FF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR, title = "lncRNA(+)", tag = c('A','B','C'))
p6_log <- create_distribution_log_plots(zscores_all %>% 
                                  dplyr::filter(Dataset %in% c("lncrna-exon1-negative-control", "lncrna-exon2-negative-control")) %>% 
                                  dplyr::select(all_of(PLOT_SELECTED_LOG)),
                                "zscores_lncRNA_negative", color = "#5271ffFF", 
                                output_dir = DISTRIBUTION_PLOTS_DIR,
                                type = "Negative Control", title = "lncRNA(-)", tag = c('D','E','F'))




# 2. Define the output directory
output_dir <- DISTRIBUTION_PLOTS_DIR

# Ensure directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save a unified plot to both PDF and PNG (same name, both extensions).
save_both <- function(plot, output_dir, name) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(output_dir, ext, paste0(name, ".", ext)),
           plot = plot, scale = 3, width = 10, height = 6, bg = "white", dpi = 300)
  }
}

# 3. Create the Function
generate_and_save_plots <- function(feature_name, display_title) {
  
  message(paste("Processing:", feature_name, "| Title:", display_title))
  
  if (is.null(p1[[feature_name]])) {
    warning(paste("Feature", feature_name, "not found in data. Skipping."))
    return(NULL)
  }
  
  # Use 'display_title' for the plot header
  common_annotation <- plot_annotation(
    title = display_title,
    theme = theme(
      plot.title = element_text(size = 84, hjust = 0.5)
    )
  )
  
  # --- Construct Plots ---
  # We still use 'feature_name' to grab the data from the lists (p1, p2, etc)
  
  plot_mrna <- (p1[[feature_name]] / p2[[feature_name]]) + common_annotation
  plot_snc  <- (p3[[feature_name]] / p4[[feature_name]]) + common_annotation
  plot_lnc  <- (p5[[feature_name]] / p6[[feature_name]]) + common_annotation
  
  # --- Save Plots ---
  # We use 'feature_name' for the file path to avoid spaces/special chars in filenames
  
  save_both(plot_mrna, output_dir, paste0("mRNA_", feature_name))
  save_both(plot_snc,  output_dir, paste0("sncRNA_", feature_name))
  save_both(plot_lnc,  output_dir, paste0("lncRNA_", feature_name))
}

generate_and_save_plots_log <- function(feature_name, display_title) {
  
  message(paste("Processing:", feature_name, "| Title:", display_title))
  
  if (is.null(p1_log[[feature_name]])) {
    warning(paste("Feature", feature_name, "not found in data. Skipping."))
    return(NULL)
  }
  
  # Use 'display_title' for the plot header
  common_annotation <- plot_annotation(
    title = display_title,
    theme = theme(
      plot.title = element_text(size = 84, hjust = 0.5)
    )
  )
  
  # --- Construct Plots ---
  # We still use 'feature_name' to grab the data from the lists (p1, p2, etc)
  
  plot_mrna <- (p1_log[[feature_name]] / p2_log[[feature_name]]) + common_annotation
  plot_snc  <- (p3_log[[feature_name]] / p4_log[[feature_name]]) + common_annotation
  plot_lnc  <- (p5_log[[feature_name]] / p6_log[[feature_name]]) + common_annotation
  
  # --- Save Plots ---
  # We use 'feature_name' for the file path to avoid spaces/special chars in filenames
  
  save_both(plot_mrna, output_dir, paste0("mRNA_", feature_name, "_log"))
  save_both(plot_snc,  output_dir, paste0("sncRNA_", feature_name, "_log"))
  save_both(plot_lnc,  output_dir, paste0("lncRNA_", feature_name, "_log"))
}

# 4. Execute the loop
# This applies the function to every item in your list
#lapply(PLOT_SELECTED_FEATURES, generate_and_save_plots)

#lapply(PLOT_SELECTED_LOG, generate_and_save_plots)

# Loop through the indices
for (i in seq_along(PLOT_SELECTED_FEATURES)) {
  
  current_feature <- PLOT_SELECTED_FEATURES[i]
  current_label   <- PLOT_SELECTED_FEATURES_LABELS[i]
  
  generate_and_save_plots(current_feature, current_label)
}

# Loop through the indices
for (i in seq_along(PLOT_SELECTED_LOG)) {
  
  current_feature <- PLOT_SELECTED_LOG[i]
  current_label   <- PLOT_SELECTED_LOG_LABELS[i]
  
  generate_and_save_plots_log(current_feature, current_label)
}


message("All plots generated.")

# Copy selected plots to Paper figures with figure-numbered names
for (ext in c("png", "pdf")) {
  file.copy(
    from      = file.path(DISTRIBUTION_PLOTS_DIR, ext, paste0("sncRNA_RPKM_primary.cell_log.", ext)),
    to        = file.path(VIOLIN_PLOTS_DIR, ext, paste0("figureS8-sncRNA-primary-cell-RPKM-z-score.", ext)),
    overwrite = TRUE
  )
  file.copy(
    from      = file.path(DISTRIBUTION_PLOTS_DIR, ext, paste0("mRNA_phyloP_max_241w.", ext)),
    to        = file.path(VIOLIN_PLOTS_DIR, ext, paste0("figureS9-mRNA-PhyloP-mammals-z-score.", ext)),
    overwrite = TRUE
  )
  file.copy(
    from      = file.path(DISTRIBUTION_PLOTS_DIR, ext, paste0("mRNA_methylome.", ext)),
    to        = file.path(VIOLIN_PLOTS_DIR, ext, paste0("figureS10-mRNA-methylome-z-score.", ext)),
    overwrite = TRUE
  )
  file.copy(
    from      = file.path(DISTRIBUTION_PLOTS_DIR, ext, paste0("mRNA_copy_number_log.", ext)),
    to        = file.path(VIOLIN_PLOTS_DIR, ext, paste0("figureS11-mRNA-copies-z-score.", ext)),
    overwrite = TRUE
  )
  file.copy(
    from      = file.path(DISTRIBUTION_PLOTS_DIR, ext, paste0("lncRNA_copy_number_log.", ext)),
    to        = file.path(VIOLIN_PLOTS_DIR, ext, paste0("figureS12-lncRNA-copies-z-score.", ext)),
    overwrite = TRUE
  )
}
