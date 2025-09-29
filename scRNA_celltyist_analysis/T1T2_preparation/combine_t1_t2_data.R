# T1/T2 Combined Data Preparation for DIVAS Analysis
# Extracts 6-block data (4 cell types + 2 bulk omics) for 114 samples

library(dplyr)

cat("=== T1/T2 Combined Analysis ===\n")

setwd("/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/scRNA_celltyist_analysis")

combined_data_dir <- "DIVAS_run/divas_input_T1T2combined"
if (dir.exists(combined_data_dir) && length(list.files(combined_data_dir, pattern = "Combined_.*\\.csv")) > 0) {
  cat("Final combined data already exists. Use: Rscript DIVAS_run/run_combined_6block_divas.R\n")
  quit(save = "no")
}
dir.create(combined_data_dir, recursive = TRUE, showWarnings = FALSE)

annotated_dir <- "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_Analysis_Repo/sc_celltype_celltypist/celltypist_server/Adult_COVID19_PBMC/annotated_data_majority_voting"
mapping_file <- "../preprocessing/process_bulk/sample_ids.tsv"
proteomics_file <- "../preprocessing/processed_omics_120/proteomics_120patients.csv"
metabolomics_file <- "../preprocessing/processed_omics_120/metabolomics_120patients.csv"

cat("Loading sample selection...\n")

# Read precomputed balanced selection
sel_path <- "T1T2_preparation/balanced_sample_selection_114.csv"
if (!file.exists(sel_path)) sel_path <- "balanced_sample_selection_114.csv"
if (!file.exists(sel_path)) stop("Missing balanced_sample_selection_114.csv")
sample_selection <- read.csv(sel_path, stringsAsFactors = FALSE)

if (!file.exists(mapping_file)) stop("Missing sample_ids.tsv file: ", mapping_file)
id_mapping <- read.csv(mapping_file, sep = "\t", stringsAsFactors = FALSE)

sample_to_incov <- setNames(id_mapping$original_filename, id_mapping$sample_id)
incov_to_sample <- setNames(id_mapping$sample_id, id_mapping$original_filename)

t1_samples <- sample_selection[sample_selection$timepoint == "T1", "sample_id"]
t2_samples <- sample_selection[sample_selection$timepoint == "T2", "sample_id"]

data_blocks <- c("CD4_T_combined", "CD8_T_combined", "CD14_Monocyte", "NK", "proteomics", "metabolomics")

celltype_map <- list(
  CD4_T_combined = "^CD4",
  CD8_T_combined = "^CD8", 
  CD14_Monocyte = "CD14 Monocyte",
  NK = "NK"
)

cat("Finding common genes with quality filtering...\n")

# Only compute common genes once for all scRNA blocks
if (!exists("global_common_genes")) {
  py_script_genes <- sprintf('
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
import math

# Load sample mapping
mapping_file = "%s"
id_mapping = pd.read_csv(mapping_file, sep="\\t")
sample_to_incov = dict(zip(id_mapping["sample_id"], id_mapping["original_filename"]))

# Target samples
all_samples = %s
annotated_dir = "%s"

# First pass: find intersection of all genes
common_genes = None
sample_gene_data = {}

print("Step 1: Finding gene intersection across all samples...")

for sample_id in all_samples:
    if sample_id not in sample_to_incov:
        continue
        
    incov_filename = sample_to_incov[sample_id]
    h5ad_file = os.path.join(annotated_dir, f"{incov_filename}_annotated.h5ad")
    
    if not os.path.exists(h5ad_file):
        continue
        
    try:
        adata = sc.read_h5ad(h5ad_file, backed="r")
        current_genes = set(adata.var_names)
        
        if common_genes is None:
            common_genes = current_genes
        else:
            common_genes = common_genes.intersection(current_genes)
            
        # Store for quality filtering step
        sample_gene_data[sample_id] = {
            "h5ad_file": h5ad_file,
            "genes": current_genes
        }
        
    except Exception as e:
        print(f"Error reading {h5ad_file}: {e}", file=sys.stderr)
        continue

if common_genes is None or len(common_genes) == 0:
    print("No common genes found across samples")
    sys.exit(1)

common_genes_list = sorted(list(common_genes))
print(f"Found {len(common_genes_list)} genes common across all samples")

# Step 2: Apply 10%% expression threshold filtering (matching original analysis)
print("Step 2: Applying 10%% expression threshold filtering...")

min_expression_percentage = 0.10  # 10%% threshold like original analysis
min_samples_for_expression = math.ceil(len(all_samples) * min_expression_percentage)

print(f"Gene must be expressed in at least {min_samples_for_expression} samples to be retained")

# Create aggregated expression matrix for quality filtering
gene_expression_counts = {}

for sample_id in all_samples:
    if sample_id not in sample_gene_data:
        continue
        
    h5ad_file = sample_gene_data[sample_id]["h5ad_file"]
    
    try:
        adata = sc.read_h5ad(h5ad_file)
        adata_common = adata[:, common_genes_list].copy()
        
        # Calculate mean expression per gene for this sample
        if adata_common.X.shape[0] > 0:  # Check if there are cells
            mean_expr = np.array(adata_common.X.mean(axis=0)).flatten()
            
            # Count genes with non-zero expression
            for i, gene in enumerate(common_genes_list):
                if gene not in gene_expression_counts:
                    gene_expression_counts[gene] = 0
                if mean_expr[i] > 0:
                    gene_expression_counts[gene] += 1
                    
    except Exception as e:
        print(f"Error processing {h5ad_file} for quality filtering: {e}", file=sys.stderr)
        continue

# Filter genes based on expression threshold
quality_filtered_genes = []
for gene in common_genes_list:
    if gene in gene_expression_counts and gene_expression_counts[gene] >= min_samples_for_expression:
        quality_filtered_genes.append(gene)

print(f"After 10%% expression filtering: {len(quality_filtered_genes)} genes retained")
print(f"Removed {len(common_genes_list) - len(quality_filtered_genes)} genes ({((len(common_genes_list) - len(quality_filtered_genes))/len(common_genes_list)*100):.1f}%%)")

# Save to file
with open("common_genes.txt", "w") as f:
    for gene in quality_filtered_genes:
        f.write(gene + "\\n")
', 
    mapping_file,
    paste0("['", paste(c(t1_samples, t2_samples), collapse = "', '"), "']"),
    annotated_dir)
  
  writeLines(py_script_genes, "find_common_genes.py")
  system("python3 find_common_genes.py")
  
  # Read quality-filtered common genes
  if (file.exists("common_genes.txt")) {
    global_common_genes <- readLines("common_genes.txt")
    file.remove(c("find_common_genes.py", "common_genes.txt"))
  } else {
    stop("Failed to find quality-filtered common genes")
  }
  
  cat("Found", length(global_common_genes), "quality-filtered common genes (10% expression threshold)\n")
}

# --- 5. Function to Extract and Combine Data ---
extract_and_combine_block <- function(block_name, t1_samples, t2_samples) {
  
  if (block_name %in% names(celltype_map)) {
    # scRNA cell-type blocks: read directly from .h5ad files
    celltype_label <- celltype_map[[block_name]]
    
    # Use Python to extract cell type data from .h5ad files using pre-computed common genes
    py_script <- sprintf('
import scanpy as sc
import pandas as pd
import os
import sys

# Load sample mapping
mapping_file = "%s"
id_mapping = pd.read_csv(mapping_file, sep="\\t")
sample_to_incov = dict(zip(id_mapping["sample_id"], id_mapping["original_filename"]))

# Target samples
t1_samples = %s
t2_samples = %s
celltype_label = "%s"
annotated_dir = "%s"

# Load pre-computed common genes
common_genes = %s

def extract_celltype_data(samples, timepoint):
    """Extract cell type data using pre-computed common genes"""
    all_data = []
    sample_names = []
    
    print(f"Extracting {celltype_label} data for {len(samples)} {timepoint} samples using {len(common_genes)} common genes...")
    
    for sample_id in samples:
        if sample_id not in sample_to_incov:
            print(f"Warning: {sample_id} not found in mapping", file=sys.stderr)
            continue
            
        incov_filename = sample_to_incov[sample_id]  # Get INCOVxxx-BL or INCOVxxx-AC
        h5ad_file = os.path.join(annotated_dir, f"{incov_filename}_annotated.h5ad")
        
        if not os.path.exists(h5ad_file):
            print(f"Warning: {h5ad_file} not found", file=sys.stderr)
            continue
            
        try:
            adata = sc.read_h5ad(h5ad_file)
            
            # Subset to common genes first (like original script)
            adata_common = adata[:, common_genes].copy()
            
            # Filter for target cell type (exact match or pattern)
            if "majority_voting" in adata_common.obs.columns:
                if celltype_label.startswith("^"):
                    # Use pattern matching for combined types
                    celltype_mask = adata_common.obs["majority_voting"].str.contains(celltype_label, na=False)
                else:
                    # Use exact matching for single types
                    celltype_mask = adata_common.obs["majority_voting"] == celltype_label
            elif "celltypist_majority_voting" in adata_common.obs.columns:
                if celltype_label.startswith("^"):
                    celltype_mask = adata_common.obs["celltypist_majority_voting"].str.contains(celltype_label, na=False)
                else:
                    celltype_mask = adata_common.obs["celltypist_majority_voting"] == celltype_label
            else:
                print(f"Warning: No cell type annotation found in {h5ad_file}", file=sys.stderr)
                continue
            
            if celltype_mask.sum() == 0:
                print(f"Warning: No {celltype_label} cells found in {sample_id}", file=sys.stderr)
                continue
            
            # Extract and aggregate (mean) expression for this cell type using common genes
            mean_expr = adata_common[celltype_mask, :].X.mean(axis=0)
            if hasattr(mean_expr, "A1"):  # Handle sparse matrices
                mean_expr = mean_expr.A1
            
            all_data.append(mean_expr)
            sample_names.append(sample_id)
            
        except Exception as e:
            print(f"Error processing {h5ad_file}: {e}", file=sys.stderr)
            continue
    
    if len(all_data) == 0:
        return pd.DataFrame()
    
    # Create DataFrame with common genes
    expr_df = pd.DataFrame(all_data).T
    expr_df.columns = sample_names
    expr_df.index = common_genes
    
    return expr_df

# Extract T1 and T2 data using pre-computed common genes
t1_data = extract_celltype_data(t1_samples, "T1")
t2_data = extract_celltype_data(t2_samples, "T2")

print(f"T1 data shape: {t1_data.shape}")
print(f"T2 data shape: {t2_data.shape}")

# Save temporary files
t1_data.to_csv("temp_t1.csv")
t2_data.to_csv("temp_t2.csv")
', 
      mapping_file,
      paste0("['", paste(t1_samples, collapse = "', '"), "']"),
      paste0("['", paste(t2_samples, collapse = "', '"), "']"),
      celltype_label,
      annotated_dir,
      paste0("['", paste(global_common_genes, collapse = "', '"), "']"))
    
    # Write and run Python script
    writeLines(py_script, "temp_extract.py")
    system("python3 temp_extract.py")
    
    # Read results
    if (file.exists("temp_t1.csv") && file.exists("temp_t2.csv")) {
      t1_data <- read.csv("temp_t1.csv", row.names = 1, check.names = FALSE)
      t2_data <- read.csv("temp_t2.csv", row.names = 1, check.names = FALSE)
      
      # Clean up temp files
      file.remove(c("temp_extract.py", "temp_t1.csv", "temp_t2.csv"))
    } else {
      stop("Failed to extract ", celltype_label, " data")
    }
    
  } else if (block_name == "proteomics") {
    if (!file.exists(proteomics_file)) stop("Proteomics file not found: ", proteomics_file)
    df <- read.csv(proteomics_file, row.names = 1, check.names = FALSE)
    t1_data <- df[, intersect(t1_samples, colnames(df)), drop = FALSE]
    t2_data <- df[, intersect(t2_samples, colnames(df)), drop = FALSE]
  } else if (block_name == "metabolomics") {
    if (!file.exists(metabolomics_file)) stop("Metabolomics file not found: ", metabolomics_file)
    df <- read.csv(metabolomics_file, row.names = 1, check.names = FALSE)
    t1_data <- df[, intersect(t1_samples, colnames(df)), drop = FALSE]
    t2_data <- df[, intersect(t2_samples, colnames(df)), drop = FALSE]
  } else {
    stop("Unknown block name: ", block_name)
  }
  
  cat("    - T1 data dimensions:", dim(t1_data), "\n")
  cat("    - T2 data dimensions:", dim(t2_data), "\n")
  
  # Check if samples exist in the data
  t1_samples_available <- intersect(t1_samples, colnames(t1_data))
  t2_samples_available <- intersect(t2_samples, colnames(t2_data))
  
  
  # Extract selected samples
  t1_extracted <- t1_data[, t1_samples_available, drop = FALSE]
  t2_extracted <- t2_data[, t2_samples_available, drop = FALSE]
  
  common_features <- intersect(rownames(t1_extracted), rownames(t2_extracted))
  
  if (length(common_features) == 0) {
    stop("No common features found between T1 and T2 for ", block_name)
  }
  
  t1_common <- t1_extracted[common_features, , drop = FALSE]
  t2_common <- t2_extracted[common_features, , drop = FALSE]
  
  combined_data <- cbind(t1_common, t2_common)
  
  return(combined_data)
}

cat("Processing data blocks...\n")

combined_data_list <- list()

for (block_name in data_blocks) {
  combined_block <- extract_and_combine_block(block_name, t1_samples, t2_samples)
  combined_data_list[[block_name]] <- combined_block
  
  output_file <- file.path(combined_data_dir, paste0("Combined_", block_name, ".csv"))
  write.csv(combined_block, output_file)
}

target_order <- c(t1_samples, t2_samples)
missing_report <- list()

for (block_name in names(combined_data_list)) {
  present <- colnames(combined_data_list[[block_name]])
  keep <- target_order[target_order %in% present]
  missing <- setdiff(target_order, keep)
  missing_report[[block_name]] <- missing
  
  combined_data_list[[block_name]] <- combined_data_list[[block_name]][, keep, drop = FALSE]
  
  output_file <- file.path(combined_data_dir, paste0("Combined_", block_name, ".csv"))
  write.csv(combined_data_list[[block_name]], output_file)
}

# Write missing report
miss_path <- file.path(combined_data_dir, "missing_samples_report.csv")
sink(miss_path)
cat("block,missing_sample\n")
for (bn in names(missing_report)) {
  if (length(missing_report[[bn]]) > 0) {
    for (ms in missing_report[[bn]]) cat(sprintf("%s,%s\n", bn, ms))
  }
}
sink()

########################################
# 7. Create Final Sample Information
########################################
final_sample_info <- sample_selection
final_sample_info$divas_sample_id <- final_sample_info$sample_id
write.csv(final_sample_info, file.path(combined_data_dir, "final_sample_info.csv"), row.names = FALSE)

data_summary <- data.frame(
  Data_Block = names(combined_data_list),
  Features = sapply(combined_data_list, nrow),
  Samples = sapply(combined_data_list, ncol),
  T1_Samples = sapply(combined_data_list, function(x) sum(colnames(x) %in% t1_samples)),
  T2_Samples = sapply(combined_data_list, function(x) sum(colnames(x) %in% t2_samples)),
  stringsAsFactors = FALSE
)

write.csv(data_summary, file.path(combined_data_dir, "data_summary.csv"), row.names = FALSE)

na_summary <- sapply(combined_data_list, function(x) sum(is.na(x)))
if (any(na_summary > 0)) {
  cat("Warning: NA values found in some blocks\n")
}

sample_orders_consistent <- all(sapply(combined_data_list[-1], function(x) {
  identical(colnames(x), colnames(combined_data_list[[1]]))
}))

# Create a README file
readme_content <- paste0(
  "# Combined T1/T2 6-Block DIVAS Analysis Data\n\n",
  "## Overview\n",
  "This directory contains the combined T1/T2 dataset for 6-block DIVAS analysis.\n",
  "Total samples: ", nrow(final_sample_info), " (balanced by WHO Ordinal Scale)\n",
  "T1 samples: ", sum(final_sample_info$timepoint == "T1"), "\n",
  "T2 samples: ", sum(final_sample_info$timepoint == "T2"), "\n\n",
  "## Data Blocks\n",
  paste(paste0("- ", data_summary$Data_Block, ": ", data_summary$Features, " features x ", data_summary$Samples, " samples"), collapse = "\n"),
  "\n\n## Files\n",
  "- Combined_*.csv: Individual data blocks\n",
  "- final_sample_info.csv: Sample metadata with WOS, demographics\n",
  "- data_summary.csv: Summary of all data blocks\n",
  "- README.md: This file\n\n",
  "Generated on: ", Sys.time(), "\n"
)

writeLines(readme_content, file.path(combined_data_dir, "README.md"))

cat("Data extraction complete -", nrow(final_sample_info), "samples\n")
cat("Next: Rscript DIVAS_run/run_combined_6block_divas.R\n")
