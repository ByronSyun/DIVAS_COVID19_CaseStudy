#!/usr/bin/env python3
# Script to convert single-cell protein expression data (CITE-seq) 
# into a pseudo-bulk format with CLR normalization for DIVAS analysis.

import os
import pandas as pd
import numpy as np
import re
import glob
from scipy.stats import gmean

# --- Path Definitions ---
# It reads raw data from the 'data_acquisition' directory and writes output to '../processed_omics_all/'.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SOURCE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../../data_acquisition/arrayexpress_data/pro_data_unzipped"))
OUTPUT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../processed_omics_all"))

# --- Function Definitions ---

def extract_sample_id(filename):
    """
    Extracts and standardizes a sample ID from a raw data filename.
    e.g., '.../100_1_pro_library.txt' -> 'COVID_100_T1'
    """
    base_name = os.path.basename(filename).replace(".txt", "")
    match = re.search(r'library_([^\.]+)', base_name)
    if not match:
        return base_name

    sample_id_str = match.group(1)

    if re.match(r'^\d+_[12]$', sample_id_str):
        patient_id, time_point = sample_id_str.split('_')
        return f"COVID_{patient_id}_T{time_point}"
    elif "BW" in sample_id_str or sample_id_str.startswith("BP") or sample_id_str in ["CL2", "Mix_donor1"]:
        return f"Healthy_{sample_id_str}"
    else:
        return f"Other_{sample_id_str}"

def clean_protein_name(name):
    """Cleans a protein name to make it a valid column header."""
    name = name.strip()
    name = re.sub(r'[^a-zA-Z0-9_\-]', '_', name)
    name = re.sub(r'_+', '_', name).strip('_')
    return name if name else 'unknown_protein'

def clr_transform(data_df):
    """
    Apply Centered Log-Ratio (CLR) transformation: clr(x) = log(x / g(x))
    """
    # Add pseudo-count to avoid log(0) issues
    pseudo_count = np.min(data_df[data_df > 0].min()) / 100
    if pd.isna(pseudo_count):
        pseudo_count = 1e-12
        
    data_with_pseudo = data_df + pseudo_count
    
    # Calculate geometric mean for each sample (column-wise)
    geometric_means = gmean(data_with_pseudo, axis=0)
    
    # Apply CLR transformation: log(x_i / g(x))
    clr_data = np.log(data_with_pseudo.div(geometric_means, axis=1))
    
    return clr_data

def convert_to_pseudo_bulk(file_path):
    """
    Reads a single-cell data file, calculates the mean expression for each protein
    across all cells, and returns it as a pandas Series.
    """
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        protein_names = [clean_protein_name(p) for p in df.columns]
        
        # Ensure data is numeric, coercing errors
        numeric_df = df.apply(pd.to_numeric, errors='coerce')
        
        # Calculate mean, ignoring NaNs that may have been introduced
        mean_expression = numeric_df.mean(axis=0, skipna=True)
        
        # Create a Series with cleaned protein names
        pseudo_bulk = pd.Series(mean_expression.values, index=protein_names)
        return pseudo_bulk

    except Exception as e:
        print(f"Warning: Could not process file {os.path.basename(file_path)}. Error: {e}")
        return None

# --- Main Execution ---

if __name__ == "__main__":
    print("--- Starting CITE-seq to Pseudo-Bulk Conversion ---")
    print(f"Reading raw data from: {SOURCE_DIR}")
    print(f"Writing output to: {OUTPUT_DIR}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    txt_files = sorted(glob.glob(os.path.join(SOURCE_DIR, "*.txt")))
    if not txt_files:
        print(f"Error: No .txt files found in the source directory.")
        print("Please ensure the 'data_acquisition' workflow has been run successfully.")
        exit(1)

    print(f"Found {len(txt_files)} CITE-seq files to process.")

    results = {}
    sample_mapping = {}

    for file_path in txt_files:
        sample_id = extract_sample_id(file_path)
        pseudo_bulk_series = convert_to_pseudo_bulk(file_path)
        if pseudo_bulk_series is not None:
            results[sample_id] = pseudo_bulk_series
            sample_mapping[sample_id] = os.path.basename(file_path)

    if not results:
        print("Error: No files were successfully processed. Aborting.")
        exit(1)

    # Combine all series into a single DataFrame (proteins as rows, samples as columns)
    datablock_df = pd.DataFrame(results)
    
    # Apply CLR normalization
    print("Applying CLR normalization...")
    datablock_df = clr_transform(datablock_df)

    # Create and save the sample ID mapping file
    sample_ids_df = pd.DataFrame.from_dict(sample_mapping, orient='index', columns=['original_filename'])
    sample_ids_df.index.name = 'sample_id'
    sample_ids_df.reset_index(inplace=True)
    sample_ids_df['sample_type'] = sample_ids_df['sample_id'].apply(lambda x: x.split('_')[0])
    
    sample_ids_path = os.path.join(OUTPUT_DIR, "sample_ids.tsv")
    sample_ids_df.to_csv(sample_ids_path, sep='\t', index=False)

    # Save the main datablock file (CLR-normalized)
    datablock_path = os.path.join(OUTPUT_DIR, "datablock_pro.tsv")
    datablock_df.to_csv(datablock_path, sep='\t')

    print("\n--- Processing Complete ---")
    print(f"Created CLR-normalized datablock with {datablock_df.shape[1]} samples and {datablock_df.shape[0]} proteins.")
    print(f"Output datablock saved to: {datablock_path}")
    print(f"Sample ID mapping file saved to: {sample_ids_path}")

    # Print a summary of sample types
    print("\nSample Type Distribution:")
    print(sample_ids_df['sample_type'].value_counts()) 