# run_combined_6block_divas.R
# DIVAS analysis on combined T1/T2 6-block dataset

library(devtools)
library(DIVAS)
library(CVXR)

cat("=== DIVAS Combined T1/T2 Analysis ===\n")

# Set working directory
setwd("/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/scRNA_celltyist_analysis")

# Define paths
input_dir <- "DIVAS_run/divas_input_T1T2combined"
output_dir <- "DIVAS_run/divas_results"
output_fig_dir <- file.path(output_dir, "diagnostic_plots")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
data_blocks <- c("CD4_T_combined", "CD8_T_combined", "CD14_Monocyte", "NK", "proteomics", "metabolomics")
data_files <- paste0("Combined_", data_blocks, ".csv")
file_paths <- file.path(input_dir, data_files)

if (!all(file.exists(file_paths))) {
  missing_files <- data_files[!file.exists(file_paths)]
  stop("Missing data files: ", paste(missing_files, collapse = ", "))
}

cat("Loading data blocks...\n")
data_frames_list <- lapply(file_paths, function(f) {
  read.csv(f, row.names = 1, check.names = FALSE)
})
names(data_frames_list) <- data_blocks

# Remove NA values and convert to matrices
data_frames_clean <- lapply(data_frames_list, function(df) {
  na_rows <- apply(df, 1, function(x) any(is.na(x)))
  if (sum(na_rows) > 0) {
    df <- df[!na_rows, , drop = FALSE]
  }
  return(df)
})

data_list <- lapply(data_frames_clean, as.matrix)

# Verify sample consistency
sample_names_per_block <- lapply(data_list, colnames)
common_samples_final <- Reduce(intersect, sample_names_per_block)

if (length(common_samples_final) != ncol(data_list[[1]])) {
  data_list <- lapply(data_list, function(x) x[, common_samples_final, drop = FALSE])
}

cat("Final data dimensions:\n")
for (name in names(data_list)) {
  cat("  ", name, ":", paste(dim(data_list[[name]]), collapse = " x "), "\n")
}

# Run DIVAS analysis
cat("Running DIVAS analysis...\n")
start_time <- Sys.time()

set.seed(123)

divas_results <- DIVASmain(
  datablock = data_list,
  nsim = 500, 
  iprint = TRUE,
  colCent = FALSE,  
  rowCent = TRUE,
  seed = 123,
  ReturnDetail = TRUE  
)

end_time <- Sys.time()
analysis_time <- end_time - start_time
cat("DIVAS analysis completed in", round(analysis_time, 2), attr(analysis_time, "units"), "\n")

# Save results
results_file <- file.path(output_dir, "divas_results_combined_6block.rds")
saveRDS(divas_results, results_file)
cat("Results saved to:", results_file, "\n")

# Generate diagnostic plots
tryCatch({
  pdf_file <- file.path(output_fig_dir, "divas_diagnostic_plots.pdf")
  pdf(pdf_file, width = 12, height = 8)
  
  DJIVEAngleDiagnosticJP(
    datablock = data_list,
    dataname = names(data_list),
    outstruct = divas_results,
    randseed = 123,
    titlestr = "DIVAS: Combined T1/T2 6-Block Analysis"
  )
  
  dev.off()
  cat("Diagnostic plots saved to:", pdf_file, "\n")
  
}, error = function(e) {
  cat("Warning: Diagnostic plot generation failed:", as.character(e), "\n")
})

cat("=== Analysis Complete ===\n")