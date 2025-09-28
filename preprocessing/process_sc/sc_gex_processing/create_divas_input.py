#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integrate single-cell GEX data for DIVAS input format
Combines processed single-cell data into gene x sample matrix for DIVAS analysis
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import glob
from datetime import datetime
import json
import sys

# Path configuration
WORKSPACE_ROOT = "/Users/byronsun/Desktop/DIVAS-code"
BASE_DIR = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing/process_sc/sc_gex_processing")
PROCESSED_DATA_DIR = os.path.join(BASE_DIR, "processed_gex_data")
OUTPUT_DIR = BASE_DIR

# Sample ID mapping file
SAMPLE_MAP_FILE = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing/process_bulk/sample_ids.tsv")

os.makedirs(OUTPUT_DIR, exist_ok=True)
log_file = os.path.join(OUTPUT_DIR, "create_divas_input.log")

def log_message(message):
    """Log message to file and console"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)

def create_reverse_naming_map(map_file_path):
    """Create reverse mapping from INCOV (original_filename) to COVID (sample_id)"""
    try:
        log_message(f"Creating reverse ID mapping from '{os.path.basename(map_file_path)}'...")
        mapping_df = pd.read_csv(map_file_path, sep='\t')
        covid_samples = mapping_df[mapping_df['sample_type'] == 'COVID'].copy()
        # Key: original_filename (INCOVxxx-BL/AC), Value: sample_id (COVID_x_Tx)
        id_map = pd.Series(covid_samples.sample_id.values, index=covid_samples.original_filename).to_dict()
        log_message(f"Successfully created ID mapping for {len(id_map)} samples.")
        return id_map
    except Exception as e:
        log_message(f"Error: Failed to create ID mapping - {e}")
        return None

def get_sample_list():
    """Get list of processed samples by scanning directory (still in INCOV format)"""
    if not os.path.isdir(PROCESSED_DATA_DIR):
        log_message(f"Error: Processed data directory '{PROCESSED_DATA_DIR}' does not exist.")
        return []
    
    # Sample IDs are subdirectory names
    sample_ids = [d for d in os.listdir(PROCESSED_DATA_DIR) if os.path.isdir(os.path.join(PROCESSED_DATA_DIR, d))]
    # Sort for deterministic processing order
    return sorted(sample_ids)

def aggregate_single_cell_data(sample_list, id_map, aggregation_method="mean"):
    """
    Aggregate single-cell data into gene x sample matrix with unified IDs
    
    Parameters:
    - sample_list: List of sample IDs (INCOV format)
    - id_map: Dictionary mapping INCOV to COVID format
    - aggregation_method: "mean" for average, "sum" for total
    
    Returns:
    - Gene x sample matrix (with unified COVID IDs)
    - Sample metadata (with unified COVID IDs)
    """
    log_message(f"Starting single-cell data aggregation using {aggregation_method} method")

    # Memory optimization strategy
    # Phase 1: Scan all samples to build gene union
    log_message("Phase 1: Scanning all samples to build gene union...")
    all_genes_set = set()
    valid_sample_paths = {}
    for sample_id in sample_list:
        sample_dir = os.path.join(PROCESSED_DATA_DIR, sample_id)
        h5ad_file = os.path.join(sample_dir, f"{sample_id}_processed.h5ad")
        if os.path.exists(h5ad_file):
            try:
                adata = sc.read_h5ad(h5ad_file, backed='r')  # Read-only mode for faster access
                all_genes_set.update(adata.var_names.tolist())
                valid_sample_paths[sample_id] = h5ad_file
            except Exception as e:
                log_message(f"Warning: Error reading gene list from sample {sample_id}: {e}")
        else:
            log_message(f"Warning: h5ad file not found for sample {sample_id}: {h5ad_file}")
            
    all_genes_sorted = sorted(list(all_genes_set))
    log_message(f"Gene union built: {len(all_genes_sorted)} unique genes found.")

    # Phase 2: Process each sample and fill final matrix
    log_message("Phase 2: Processing samples and filling final data matrix...")
    
    # Initialize empty DataFrame with full gene set as index
    aggregated_matrix = pd.DataFrame(index=all_genes_sorted)
    sample_metadata = {}

    for sample_id, h5ad_file in valid_sample_paths.items():  # sample_id still in INCOV format
        try:
            # Core ID conversion
            unified_sample_id = id_map.get(sample_id, sample_id)
            log_message(f"Processing sample: {sample_id} -> {unified_sample_id}")

            # Check for duplicate IDs to prevent data overwrite
            if unified_sample_id in aggregated_matrix.columns:
                log_message(f"Warning: Sample ID '{unified_sample_id}' is duplicate. Source file {sample_id} will be skipped.")
                continue

            adata = sc.read_h5ad(h5ad_file)  # Full read
            
            # Calculate pseudo-bulk expression
            if aggregation_method == "mean":
                expr_vector = adata.X.mean(axis=0)
            elif aggregation_method == "sum":
                expr_vector = adata.X.sum(axis=0)
            
            # Create Series with current sample genes as index
            sample_series = pd.Series(
                expr_vector.A1 if hasattr(expr_vector, 'A1') else expr_vector,
                index=adata.var_names
            )

            # Add Series to main matrix, pandas will auto-align by index
            # Genes in main matrix but not in current sample will be filled with NaN
            aggregated_matrix[unified_sample_id] = sample_series

            # Record metadata using unified ID as key
            condition = "COVID" if "INCOV" in sample_id else "Healthy"
            is_baseline = "-BL" in sample_id
            time_point = "Baseline" if is_baseline else "After COVID"
            sample_metadata[unified_sample_id] = {
                "condition": condition,
                "time_point": time_point,
                "cell_count": adata.n_obs
            }
            log_message(f"Successfully aggregated sample {unified_sample_id} with {adata.n_obs} cells")

        except Exception as e:
            log_message(f"Error processing sample {sample_id}: {str(e)}")

    # Fill missing values with 0
    log_message("Filling all missing values (NaN) with 0...")
    aggregated_matrix.fillna(0, inplace=True)
    
    # Create sample metadata dataframe
    metadata_df = pd.DataFrame.from_dict(sample_metadata, orient='index')
    
    return aggregated_matrix, metadata_df

def main():
    """Main function"""
    log_message("Starting DIVAS input data creation (with built-in ID unification)")
    
    # 1. Create ID mapping
    id_map = create_reverse_naming_map(SAMPLE_MAP_FILE)
    if id_map is None:
        log_message("Error: Unable to create ID mapping, terminating.")
        sys.exit(1)

    # 2. Get sample list (INCOV format)
    sample_list = get_sample_list()
    log_message(f"Found {len(sample_list)} samples to process")
    
    # 3. Aggregate single-cell data with ID conversion
    divas_datablock, metadata = aggregate_single_cell_data(sample_list, id_map, aggregation_method="mean")
    
    # Save as DIVAS format datablock (genes x samples)
    datablock_file = os.path.join(OUTPUT_DIR, "gex_divas_datablock.tsv")
    divas_datablock.to_csv(datablock_file, sep='\t')
    log_message(f"Saved DIVAS datablock to: {datablock_file}, dimensions: {divas_datablock.shape}")
    
    # Save sample metadata
    metadata_file = os.path.join(OUTPUT_DIR, "gex_sample_metadata.tsv")
    metadata.to_csv(metadata_file, sep='\t')
    log_message(f"Saved sample metadata to: {metadata_file}, dimensions: {metadata.shape}")
    
    log_message("DIVAS input data creation complete")

if __name__ == "__main__":
    main()