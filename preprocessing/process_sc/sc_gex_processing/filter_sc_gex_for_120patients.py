#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter scRNA-seq raw data directory to keep only files for 120 dual time-point patients
"""

import os
import pandas as pd
import glob
import re

def get_target_standard_ids():
    """
    Determine the 240 standardized sample IDs (e.g., COVID_1_T1) for 120 dual time-point patients
    """
    # Define base paths
    WORKSPACE_ROOT = "/Users/byronsun/Desktop/DIVAS-code"
    BASE_DIR = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing")
    
    # Read core samples and ID mapping table
    core_samples_file = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/sample_distribution/core_samples_10X_metabolomics_proteomics_20250704_173841.csv"
    mapping_file = os.path.join(BASE_DIR, "process_bulk/sample_ids.tsv")
    
    core_samples = pd.read_csv(core_samples_file)
    id_mapping = pd.read_csv(mapping_file, sep='\t')

    # Filter for 120 dual time-point patients
    patient_counts = core_samples['Study Subject ID'].value_counts()
    dual_timepoint_patients = patient_counts[patient_counts == 2].index.tolist()
    
    print(f"Found {len(dual_timepoint_patients)} dual time-point patients")
    
    # Get corresponding sample IDs
    dual_samples = core_samples[core_samples['Study Subject ID'].isin(dual_timepoint_patients)]
    target_sample_ids = dual_samples['Sample ID'].tolist()
    
    print(f"Target sample IDs: {len(target_sample_ids)}")
    
    # Convert to mapping format
    def convert_to_mapping_format(sample_id):
        if sample_id.endswith('-1'):
            return sample_id.replace('-1', '-BL')
        elif sample_id.endswith('-2'):
            return sample_id.replace('-2', '-AC')
        return sample_id
    
    mapping_format_ids = [convert_to_mapping_format(sid) for sid in target_sample_ids]
    
    # Find standardized IDs
    target_mapping = id_mapping[id_mapping['original_filename'].isin(mapping_format_ids)]
    standard_ids = target_mapping['sample_id'].tolist()
    
    print(f"Standardized sample IDs: {len(standard_ids)}")
    return standard_ids

def extract_sample_info_from_filename(filename):
    """
    Extract sample information from GEX filename
    Returns (library_id, batch) or None if parsing fails
    """
    # Pattern: heathlab_dc_9_17_pbmc_gex_library_83_1.txt
    pattern = r'heathlab_dc_9_17_pbmc_gex_library_(\w+)_(\d+)\.txt'
    match = re.match(pattern, filename)
    
    if match:
        library_id = match.group(1)
        batch = match.group(2)
        return (library_id, batch)
    
    return None

def main():
    print("Starting GEX source directory cleanup for 120 dual time-point patients...")
    
    # Get target standardized sample IDs
    target_standard_ids = get_target_standard_ids()
    
    # Source directory containing raw GEX files
    source_dir = "gex_data_unzipped"
    
    if not os.path.exists(source_dir):
        print(f"Error: Source directory '{source_dir}' not found")
        return
    
    # Get all .txt files
    all_files = glob.glob(os.path.join(source_dir, "*.txt"))
    print(f"Found {len(all_files)} total files in source directory")
    
    # Read sample ID mapping for reverse lookup
    mapping_file = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/process_bulk/sample_ids.tsv"
    id_mapping = pd.read_csv(mapping_file, sep='\t')
    
    # Create reverse mapping: standard_id -> (library_id, batch)
    standard_to_file_info = {}
    
    for _, row in id_mapping.iterrows():
        sample_id = row['sample_id']
        if sample_id in target_standard_ids:
            # Extract library_id and batch from sample_id
            if 'COVID_' in sample_id:
                parts = sample_id.split('_')
                if len(parts) >= 3:
                    try:
                        patient_num = int(parts[1])
                        time_point = parts[2]
                        batch = "1" if time_point == "T1" else "2"
                        library_id = str(patient_num)
                        standard_to_file_info[sample_id] = (library_id, batch)
                    except:
                        continue
    
    print(f"Created file info mapping for {len(standard_to_file_info)} samples")
    
    # Find files to keep
    files_to_keep = []
    target_file_patterns = set()
    
    for standard_id, (library_id, batch) in standard_to_file_info.items():
        target_file_patterns.add((library_id, batch))
    
    print(f"Looking for {len(target_file_patterns)} unique file patterns")
    
    # Check each file
    for file_path in all_files:
        filename = os.path.basename(file_path)
        file_info = extract_sample_info_from_filename(filename)
        
        if file_info and file_info in target_file_patterns:
            files_to_keep.append(file_path)
    
    print(f"Files to keep: {len(files_to_keep)}")
    print(f"Files to remove: {len(all_files) - len(files_to_keep)}")
    
    # Ask for confirmation
    response = input(f"\nProceed to remove {len(all_files) - len(files_to_keep)} files? (y/n): ")
    
    if response.lower() == 'y':
        removed_count = 0
        for file_path in all_files:
            if file_path not in files_to_keep:
                try:
                    os.remove(file_path)
                    removed_count += 1
                except Exception as e:
                    print(f"Error removing {file_path}: {e}")
        
        print(f"Successfully removed {removed_count} files")
        print(f"Remaining files: {len(files_to_keep)}")
        print("GEX source directory cleanup complete!")
    else:
        print("Cleanup cancelled")

if __name__ == "__main__":
    main()
