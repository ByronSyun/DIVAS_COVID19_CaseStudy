#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single-cell GEX data processing script
Processes and standardizes single-cell gene expression data for DIVAS analysis
"""

import os
import pandas as pd
import numpy as np
import glob
from datetime import datetime
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import re
import sys

# Path configuration
WORKSPACE_ROOT = "/Users/byronsun/Desktop/DIVAS-code"
BASE_DIR = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing/process_sc/sc_gex_processing")
OUTPUT_DIR = os.path.join(BASE_DIR, "processed_gex_data")
SAMPLE_IDS_FILE = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing/process_bulk/sample_ids.tsv")

os.makedirs(OUTPUT_DIR, exist_ok=True)
log_file = os.path.join(OUTPUT_DIR, "process_gex_data.log")

def log_message(message):
    """Log message to file and console"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)
    with open(log_file, "a") as f:
        f.write(log_msg + "\n")

def create_sample_id_mappings():
    """Create sample ID mappings from library IDs to standardized sample IDs"""
    log_message(f"Reading sample ID mapping file: {SAMPLE_IDS_FILE}")
    
    try:
        sample_ids_df = pd.read_csv(SAMPLE_IDS_FILE, sep='\t')
        
        # Create mapping dictionaries
        library_to_incov = {}  # COVID samples: (library_id, batch) -> INCOV sample ID
        healthy_ids = {}  # Healthy controls: library_id -> Healthy sample ID
        for _, row in sample_ids_df.iterrows():
            sample_id = row['sample_id']
            original_id = row['original_filename']
            sample_type = row['sample_type']
            
            # Process COVID samples (INCOV format)
            if isinstance(original_id, str) and original_id.startswith('INCOV'):
                match = re.match(r'COVID_(\d+)_T(\d+)', sample_id)
                if match:
                    num = int(match.group(1))
                    batch = int(match.group(2))
                    
                    library_id = str(num)
                    # Create mapping
                    key = (library_id, str(batch))
                    library_to_incov[key] = original_id
            
            # Process healthy control samples
            elif sample_type == 'Healthy' and sample_id.startswith('Healthy_'):
                identifier = sample_id.replace('Healthy_', '')
                healthy_ids[identifier] = sample_id
        
        log_message(f"Created {len(library_to_incov)} COVID sample ID mappings")
        log_message(f"Created {len(healthy_ids)} healthy control sample ID mappings")
        return library_to_incov, healthy_ids
    
    except Exception as e:
        log_message(f"Error creating mappings: {str(e)}")
        return {}, {}

def process_single_file(file_path, output_dir, sample_id):
    """Process single scRNA-seq data file"""
    try:
        log_message(f"Processing file: {os.path.basename(file_path)}, Sample ID: {sample_id}")
        
        # Create sample-specific output directory
        sample_output_dir = os.path.join(output_dir, sample_id)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Read single-cell data
        adata = sc.read_text(file_path, delimiter='\t')
        log_message(f"Data dimensions: {adata.shape}")
        
        # Record original sample ID
        adata.obs['sample_id'] = sample_id
        
        # Basic preprocessing
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # Calculate QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        # Visualize QC metrics
        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        sns.histplot(adata.obs.n_genes_by_counts, kde=False, ax=axs[0])
        axs[0].set_title(f'Distribution of genes per cell - {sample_id}')
        sns.histplot(adata.obs.pct_counts_mt, kde=False, ax=axs[1])
        axs[1].set_title(f'Distribution of mitochondrial gene percentage - {sample_id}')
        plt.tight_layout()
        plt.savefig(os.path.join(sample_output_dir, "qc_metrics.png"))
        plt.close()
        
        # Filter low-quality cells
        adata = adata[adata.obs.n_genes_by_counts < 2500, :]
        adata = adata[adata.obs.pct_counts_mt < 5, :]
        
        # Normalization
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Save processed data
        output_file = os.path.join(sample_output_dir, f"{sample_id}_processed.h5ad")
        adata.write(output_file)
        log_message(f"Sample {sample_id} processed and saved to {output_file}")
        
        return True, f"Sample {sample_id} processed successfully"
    
    except Exception as e:
        log_message(f"Error processing sample {sample_id}: {str(e)}")
        return False, str(e)

def main():
    """Main processing function"""
    
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
        log_message(f"Data directory from command line: {data_dir}")
    else:
        data_dir = os.path.join(BASE_DIR, "gex_data_unzipped")
        log_message(f"Using default data directory: {data_dir}")

    log_message("Starting single-cell GEX data processing")
    log_message(f"Data directory: {data_dir}")
    log_message(f"Output directory: {OUTPUT_DIR}")
    
    # Check data files
    data_files = glob.glob(os.path.join(data_dir, "*.txt"))
    log_message(f"Found {len(data_files)} data files")
    
    if len(data_files) == 0:
        log_message("Error: No data files found, please run data extraction first")
        return
    
    # Create sample ID mappings
    library_to_incov, healthy_ids = create_sample_id_mappings()
    
    # Extract sample information from filenames
    sample_info = []
    covid_count = 0
    healthy_count = 0
    
    for file in data_files:
        filename = os.path.basename(file)
        # Parse filename to get library_id and batch
        sample_id = "unknown"
        is_covid = "unknown"
        
        # Process healthy control samples
        if any(healthy_id in filename for healthy_id in healthy_ids.keys()):
            for healthy_id in healthy_ids.keys():
                if healthy_id in filename:
                    sample_id = healthy_ids[healthy_id]
                    is_covid = "Healthy"
                    library_id = healthy_id
                    healthy_count += 1
                    break
        # Process mixed samples
        elif "Mix_donor" in filename:
            parts = filename.split("_")
            donor_part = parts[-1].replace(".txt", "")
            sample_id = f"Mix_donor_{donor_part}"
            library_id = "Mix"
            is_covid = "Control"
        # Process standard COVID samples
        else:
            parts = filename.split("_")
            if len(parts) >= 9:
                batch = parts[8].replace(".txt", "")
                library_id = parts[7]
                
                # Use mapping to find INCOV sample ID
                key = (library_id, batch)
                if key in library_to_incov:
                    sample_id = library_to_incov[key]
                    is_covid = "COVID"
                    covid_count += 1
                else:
                    sample_id = f"unknown_{library_id}_{batch}"
            else:
                library_id = "unknown"
        
        sample_type = "PBMC"
        
        sample_info.append({
            "file": file,
            "filename": filename,
            "sample_id": sample_id,
            "library_id": library_id,
            "sample_type": sample_type,
            "condition": is_covid
        })
    
    sample_df = pd.DataFrame(sample_info)
    log_message(f"Extracted {len(sample_df)} sample information")
    log_message(f"COVID samples: {covid_count}, Healthy controls: {healthy_count}")
    
    # Save sample information
    sample_info_file = os.path.join(OUTPUT_DIR, "sample_info.csv")
    sample_df.to_csv(sample_info_file, index=False)
    log_message(f"Sample info saved to: {sample_info_file}")
    
    # Sort by sample ID
    log_message("Sorting by sample ID...")
    
    # Sort COVID samples by INCOV number
    covid_samples = sample_df[sample_df['condition'] == 'COVID'].copy()
    covid_samples['sort_key'] = covid_samples['sample_id'].apply(
        lambda x: int(re.search(r'INCOV(\d+)', x).group(1)) if isinstance(x, str) and re.search(r'INCOV(\d+)', x) else 999999
    )
    covid_samples = covid_samples.sort_values('sort_key')
    
    # Sort healthy control samples
    healthy_samples = sample_df[sample_df['condition'] == 'Healthy'].copy()
    
    # Other samples
    other_samples = sample_df[~sample_df['condition'].isin(['COVID', 'Healthy'])]
    
    # Combine sorted samples
    sorted_samples = pd.concat([covid_samples, healthy_samples, other_samples])
    sorted_samples = sorted_samples.drop('sort_key', axis=1, errors='ignore')
    
    # Process all found samples
    log_message(f"Starting to process {len(sorted_samples)} samples...")
    
    processed_count = 0
    for _, row in sorted_samples.iterrows():
        file_path = row['file']
        sample_id = row['sample_id']
        
        # Pass main output directory to processing function
        success, message = process_single_file(file_path, OUTPUT_DIR, sample_id)
        if success:
            processed_count += 1
    
    log_message(f"Processing complete: {processed_count}/{len(sorted_samples)} samples processed successfully")
    
    log_message("Single-cell GEX data processing finished")

if __name__ == "__main__":
    main() 