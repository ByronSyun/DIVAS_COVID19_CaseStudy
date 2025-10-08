#!/usr/bin/env python3
"""
Extract cell-level metadata from all CellTypist annotated .h5ad files.

Output CSV columns:
- cell_id, patient_id, sample_name, timepoint
- majority_voting, predicted_labels
- n_genes, total_counts, pct_counts_mt

Usage: python extract_all_cell_metadata.py
"""

import os
import scanpy as sc
import pandas as pd

def extract_all_metadata(input_dir, output_file):
    """Extract metadata from all annotated .h5ad files and combine into one CSV."""
    
    h5ad_files = sorted([f for f in os.listdir(input_dir) if f.endswith('_annotated.h5ad')])
    print(f"Found {len(h5ad_files)} annotated files")
    
    all_metadata = []
    
    for idx, h5ad_file in enumerate(h5ad_files, 1):
        sample_name = h5ad_file.replace('_annotated.h5ad', '')
        file_path = os.path.join(input_dir, h5ad_file)
        
        print(f"[{idx}/{len(h5ad_files)}] {sample_name}", end=" ")
        
        try:
            adata = sc.read_h5ad(file_path)
            df = adata.obs.copy()
            df['cell_id'] = df.index
            
            # Extract patient ID and timepoint (Format: INCOV###-BL/AC)
            parts = sample_name.split('-')
            patient_id = parts[0]
            timepoint_code = parts[1] if len(parts) > 1 else 'Unknown'
            timepoint_map = {'BL': 'T1', 'AC': 'T2'}
            timepoint = timepoint_map.get(timepoint_code, timepoint_code)
            
            df['patient_id'] = patient_id
            df['timepoint'] = timepoint
            df['sample_name'] = sample_name
            
            column_order = [
                'cell_id', 'patient_id', 'sample_name', 'timepoint',
                'majority_voting', 'predicted_labels', 'sample_id',
                'n_genes', 'n_genes_by_counts', 'total_counts',
                'total_counts_mt', 'pct_counts_mt'
            ]
            existing_columns = [col for col in column_order if col in df.columns]
            df = df[existing_columns]
            
            all_metadata.append(df)
            print(f"({len(df)} cells)")
            
        except Exception as e:
            print(f"ERROR: {e}")
            continue
    
    if all_metadata:
        print("\nCombining and saving...")
        combined_df = pd.concat(all_metadata, ignore_index=False)
        combined_df.to_csv(output_file, index=False)
        
        print(f"\nSUCCESS: {len(combined_df):,} total cells from {len(h5ad_files)} samples")
        print(f"Output: {output_file} ({os.path.getsize(output_file) / 1024 / 1024:.1f} MB)")
        print(f"\nCell type distribution:")
        print(combined_df['majority_voting'].value_counts())
        
        return combined_df
    else:
        print("\nERROR: No metadata extracted!")
        return None

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # TODO: After uploading to Zenodo, users should download and place data in:
    # annotated_data_majority_voting/
    # For now, using analysis repo data
    input_dir = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_Analysis_Repo/sc_celltype_celltypist/celltypist_server/Adult_COVID19_PBMC/annotated_data_majority_voting"
    
    # Future: Use local directory after Zenodo download
    # input_dir = os.path.join(script_dir, "annotated_data_majority_voting")
    
    if not os.path.exists(input_dir):
        print(f"ERROR: Annotated data not found at: {input_dir}")
        print(f"Please download data from Zenodo and place in annotated_data_majority_voting/")
        exit(1)
    
    output_file = os.path.join(script_dir, "all_cells_metadata_complete.csv")
    
    print("="*60)
    print("CellTypist Metadata Extraction")
    print("="*60)
    print(f"Input:  {input_dir}")
    print(f"Output: {output_file}\n")
    
    result = extract_all_metadata(input_dir, output_file)
    
    if result is not None:
        print(f"\nFirst 5 rows preview:")
        print(result.head(5).to_string())

