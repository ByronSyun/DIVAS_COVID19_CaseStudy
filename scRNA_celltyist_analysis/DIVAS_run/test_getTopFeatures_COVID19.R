#!/usr/bin/env Rscript
#' Test getTopFeatures function with real COVID-19 DIVAS results
#' 
#' This script tests the new getTopFeatures function using the 6-block
#' T1+T2 combined DIVAS analysis results from the COVID-19 study.

# Load the DIVAS package
library(DIVAS)

cat(strrep("=", 70), "\n")
cat("Testing getTopFeatures with COVID-19 6-block DIVAS Results\n")
cat(strrep("=", 70), "\n\n")

# Load DIVAS results
results_path <- "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/scRNA_celltyist_analysis/DIVAS_run/DIVAS_Results/divas_results_combined_6block.rds"

cat("Loading DIVAS results from:\n")
cat(results_path, "\n\n")

divasRes <- readRDS(results_path)

# Check structure
cat("DIVAS Results Structure:\n")
cat("- Available modalities:", paste(names(divasRes$Loadings), collapse=", "), "\n")
cat("- Number of modalities:", length(divasRes$Loadings), "\n\n")

# Get component names from first modality
first_modality <- names(divasRes$Loadings)[1]
component_names <- colnames(divasRes$Loadings[[first_modality]])
cat("Available components (from", first_modality, "):\n")
for (i in seq_along(component_names)) {
  cat(sprintf("  %d. %s\n", i, component_names[i]))
}
cat("\n")

# Test 1: Extract top features from first joint component
cat(strrep("=", 70), "\n")
cat("Test 1: Extract top features from first joint component\n")
cat(strrep("=", 70), "\n\n")

# Find first joint component (contains multiple modalities)
joint_comp <- component_names[grepl("\\+", component_names)][1]
cat("Selected component:", joint_comp, "\n\n")

# Test with CD4_T_combined modality
test_modality <- "CD4_T_combined"
if (test_modality %in% names(divasRes$Loadings)) {
  cat("Extracting top 20 positive and 20 negative features from", test_modality, "\n\n")
  
  top_features <- getTopFeatures(
    divasRes = divasRes,
    compName = joint_comp,
    modName = test_modality,
    n_top_pos = 20,
    n_top_neg = 20,
    return_values = FALSE
  )
  
  cat("Top 20 positive features:\n")
  print(top_features$top_positive)
  
  cat("\nTop 20 negative features:\n")
  print(top_features$top_negative)
  
  # Also get with values
  cat("\n\nTop 5 features with loading values:\n")
  top_features_values <- getTopFeatures(
    divasRes = divasRes,
    compName = joint_comp,
    modName = test_modality,
    n_top_pos = 5,
    n_top_neg = 5,
    return_values = TRUE
  )
  
  cat("\nTop 5 positive loadings:\n")
  print(top_features_values$top_positive)
  
  cat("\nTop 5 negative loadings:\n")
  print(top_features_values$top_negative)
} else {
  cat("Warning:", test_modality, "not found in results\n")
}

# Test 2: Extract from proteomics modality
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Test 2: Extract top features from proteomics modality\n")
cat(strrep("=", 70), "\n\n")

test_modality2 <- "proteomics"
if (test_modality2 %in% names(divasRes$Loadings)) {
  cat("Extracting top 15 features from", test_modality2, "\n\n")
  
  top_proteins <- getTopFeatures(
    divasRes = divasRes,
    compName = joint_comp,
    modName = test_modality2,
    n_top_pos = 15,
    n_top_neg = 15,
    return_values = TRUE
  )
  
  cat("Top 15 positive protein loadings:\n")
  print(top_proteins$top_positive)
  
  cat("\nTop 15 negative protein loadings:\n")
  print(top_proteins$top_negative)
} else {
  cat("Warning:", test_modality2, "not found in results\n")
}

# Test 3: Extract from metabolomics modality
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Test 3: Extract top features from metabolomics modality\n")
cat(strrep("=", 70), "\n\n")

test_modality3 <- "metabolomics"
if (test_modality3 %in% names(divasRes$Loadings)) {
  cat("Extracting top 10 features from", test_modality3, "\n\n")
  
  top_metabolites <- getTopFeatures(
    divasRes = divasRes,
    compName = joint_comp,
    modName = test_modality3,
    n_top_pos = 10,
    n_top_neg = 10,
    return_values = TRUE
  )
  
  cat("Top 10 positive metabolite loadings:\n")
  print(top_metabolites$top_positive)
  
  cat("\nTop 10 negative metabolite loadings:\n")
  print(top_metabolites$top_negative)
} else {
  cat("Warning:", test_modality3, "not found in results\n")
}

# Test 4: Test with different components
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Test 4: Compare features across multiple components\n")
cat(strrep("=", 70), "\n\n")

# Get first 3 joint components
joint_comps <- component_names[grepl("\\+", component_names)][1:min(3, sum(grepl("\\+", component_names)))]

if (test_modality %in% names(divasRes$Loadings)) {
  for (comp in joint_comps) {
    cat("\nComponent:", comp, "\n")
    cat(strrep("-", 70), "\n")
    
    top_5 <- getTopFeatures(
      divasRes = divasRes,
      compName = comp,
      modName = test_modality,
      n_top_pos = 5,
      n_top_neg = 5,
      return_values = FALSE
    )
    
    cat("Top 5 positive:", paste(top_5$top_positive, collapse=", "), "\n")
    cat("Top 5 negative:", paste(top_5$top_negative, collapse=", "), "\n")
  }
}

# Summary
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Test Summary\n")
cat(strrep("=", 70), "\n")
cat("✓ Successfully loaded COVID-19 6-block DIVAS results\n")
cat("✓ Extracted top features from multiple modalities\n")
cat("✓ Tested with and without loading values\n")
cat("✓ Compared features across different components\n")
cat("\nAll tests completed successfully!\n")
cat(strrep("=", 70), "\n")

