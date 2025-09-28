#!/usr/bin/env Rscript

# Improve Metabolomics Data Quality
#
# This script performs the following steps:
# 1. Assesses the distribution of missing values.
# 2. Filters out metabolites with a high percentage of missing values (>20%).
# 3. Imputes the remaining missing values using the K-Nearest Neighbors (KNN) method.
# 4. Applies Probabilistic Quotient Normalization (PQN) to the data.
# 5. Generates quality control plots (e.g., PCA) to evaluate the results.
# 6. Saves the cleaned, normalized data.

# Load required packages
packages <- c("impute", "pcaMethods", "pheatmap", "ggplot2", "dplyr", "tidyr")
for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        # The BiocManager installation is handled by the user or a setup script
        # if needed. This script assumes packages are available.
        stop(paste("Package", pkg, "is not installed. Please install it to continue."))
    }
}

library(impute)
library(pcaMethods)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)

cat("Starting metabolomics data quality improvement...\n")

# Define paths
# Assumes the script is run from the 'process_bulk' directory
data_dir <- "."
output_dir <- "improved_metabolomics"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data processed by the python script
data_path <- file.path(data_dir, "metabolomics_data.csv")
if (!file.exists(data_path)) {
    stop(paste("Input data file not found:", data_path))
}
metabolomics_data <- read.csv(data_path, row.names = 1, check.names = FALSE)
cat("Loaded metabolomics data. Dimensions:", dim(metabolomics_data), "\n")

# --- 1. Missing Value Assessment ---
cat("Assessing missing values...\n")
na_count_per_metabolite <- rowSums(is.na(metabolomics_data))
na_percent_per_metabolite <- na_count_per_metabolite / ncol(metabolomics_data) * 100

# --- 2. Filter Metabolites ---
missing_threshold <- 20
metabolites_to_keep <- na_percent_per_metabolite <= missing_threshold
original_metabolite_count <- nrow(metabolomics_data)
filtered_data <- metabolomics_data[metabolites_to_keep, ]
retained_count <- nrow(filtered_data)
cat("Filtered out", original_metabolite_count - retained_count, "metabolites with >", missing_threshold, "% missing values.\n")
cat("Retained", retained_count, "metabolites.\n")

# --- 3. KNN Imputation ---
cat("Imputing remaining missing values using KNN (k=10)...\n")
# impute.knn expects samples in rows and features in columns, so we transpose
transposed_data <- t(filtered_data)
knn_result <- impute.knn(transposed_data, k = 10)
imputed_data <- t(knn_result$data)
cat("KNN imputation complete. Remaining NA:", sum(is.na(imputed_data)), "\n")

if (sum(is.na(imputed_data)) > 0) {
  cat("Warning: NA values remain after imputation. Removing rows with NAs.\n")
  imputed_data <- imputed_data[complete.cases(imputed_data), ]
}

# --- 4. PQN Normalization ---
cat("Applying Probabilistic Quotient Normalization (PQN)...\n")
reference_spectrum <- apply(imputed_data, 1, median, na.rm = TRUE)
quotients <- sweep(imputed_data, 1, reference_spectrum, FUN = "/")
normalization_factors <- apply(quotients, 2, median, na.rm = TRUE)
normalized_data <- sweep(imputed_data, 2, normalization_factors, FUN = "/")
cat("PQN normalization complete.\n")

# --- 5. Save Processed Data ---
output_file_path <- file.path(output_dir, "improved_metabolomics_data.csv")
write.csv(normalized_data, output_file_path)
cat("Saved improved metabolomics data to:", output_file_path, "\n")

# --- 6. PCA for Quality Assessment ---
cat("Performing PCA for quality assessment...\n")
pca_result <- pca(t(normalized_data), nPcs = 5, scale = "uv")
score_data <- as.data.frame(pca_result@scores)
score_data$Sample <- rownames(pca_result@scores)

# Add group information for coloring
score_data$Group <- "Unknown"
score_data$Group[grepl("COVID", score_data$Sample, ignore.case = TRUE)] <- "COVID"
score_data$Group[grepl("HEALTHY|HD", score_data$Sample, ignore.case = TRUE)] <- "Healthy"

# Create PCA plot
pca_plot <- ggplot(score_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = "PCA of Improved Metabolomics Data",
       x = paste0("PC1 (", round(pca_result@R2[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(pca_result@R2[2] * 100, 1), "%)")) +
  scale_color_manual(values = c("COVID" = "darkred", "Healthy" = "darkblue", "Unknown" = "grey"))

# Save PCA plot
pca_plot_path <- file.path(output_dir, "pca_plot_metabolomics.png")
ggsave(pca_plot_path, pca_plot, width = 8, height = 6)
cat("Saved PCA plot to:", pca_plot_path, "\n")

cat("Metabolomics data quality improvement complete.\n") 