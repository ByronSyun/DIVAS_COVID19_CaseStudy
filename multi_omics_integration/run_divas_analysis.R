# DIVAS Analysis for COVID-19 4-omics Dataset
# Run DIVAS on proteomics, metabolomics, sc-GEX, and sc-protein data (120 patients)

# Set working directory to DIVAS_COVID19_CaseStudy root
setwd("/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy")
cat("Working directory:", getwd(), "\n")


# Install DIVAS from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("ByronSyun/DIVAS_Develop/pkg", ref = "main", force = TRUE)


# Load libraries
suppressPackageStartupMessages({
  library(DIVAS)
  library(ggplot2)
})

cat("Starting DIVAS analysis for 4-omics data...\n")

# Load data from processed_omics_120 directory
data_dir <- "preprocessing/processed_omics_120"

# Load data files
proteomics_df <- read.csv(file.path(data_dir, "proteomics_120patients.csv"), row.names = 1, check.names = FALSE)
metabolomics_df <- read.csv(file.path(data_dir, "metabolomics_120patients.csv"), row.names = 1, check.names = FALSE)
sc_gex_df <- read.csv(file.path(data_dir, "sc_gex_120patients_aligned.csv"), row.names = 1, check.names = FALSE)
sc_pro_df <- read.csv(file.path(data_dir, "sc_pro_120patients.csv"), row.names = 1, check.names = FALSE)

cat(sprintf("Proteomics: %d features x %d samples\n", nrow(proteomics_df), ncol(proteomics_df)))
cat(sprintf("Metabolomics: %d features x %d samples\n", nrow(metabolomics_df), ncol(metabolomics_df)))
cat(sprintf("SC-GEX: %d features x %d samples\n", nrow(sc_gex_df), ncol(sc_gex_df)))
cat(sprintf("SC-Proteomics: %d features x %d samples\n", nrow(sc_pro_df), ncol(sc_pro_df)))


# Prepare data for DIVAS
data_list <- list(
  Proteomics = as.matrix(proteomics_df),
  Metabolomics = as.matrix(metabolomics_df),
  scGEX = as.matrix(sc_gex_df),
  scProtein = as.matrix(sc_pro_df)
)


# Run DIVAS analysis with updated parameters
cat("Running DIVAS algorithm...\n")
start_time <- Sys.time()

# Create output directory  
output_dir <- "multi_omics_integration/DIVAS_Results"
output_fig_dir <- file.path(output_dir, "diagnostic_plots")
if (!dir.exists(output_fig_dir)) {
  dir.create(output_fig_dir, recursive = TRUE)
}

# Run DIVAS with latest parameters
set.seed(123)
divas_results <- DIVASmain(
  datablock = data_list,
  nsim = 400,
  iprint = TRUE,
  colCent = TRUE,
  rowCent = FALSE,
  seed = 123,
  ReturnDetail = TRUE,
  figdir = output_fig_dir
)

end_time <- Sys.time()
cat(sprintf("DIVAS analysis complete. Time elapsed: %s\n", format(end_time - start_time)))


# Save results to DIVAS_Results directory
results_path <- file.path(output_dir, "divas_results_4omics_full.rds")
saveRDS(divas_results, file = results_path)
cat(sprintf("Results saved to: %s\n", results_path))

# Generate diagnostic plots
dataname <- names(data_list)
diagnostic_plots <- DJIVEAngleDiagnosticJP(
  datablock = data_list,
  dataname = dataname,
  outstruct = divas_results,
  randseed = 566,
  titlestr = "COVID-19 4-Omics DIVAS Analysis"
)

# Save plots to PDF
plot_output_path <- file.path(output_dir, "divas_diagnostic_plots_4omics_full.pdf")
pdf(plot_output_path, width = 14, height = 9)
for (plot_item in diagnostic_plots) {
  try(print(plot_item))
}
dev.off()
cat(sprintf("Diagnostic plots saved to: %s\n", plot_output_path))

cat("DIVAS analysis completed successfully\n") 