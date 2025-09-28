#!/usr/bin/env python3
"""
Filter single-cell protein data (datablock_pro.tsv) for 120 dual time-point patients
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

def main():
    print("Starting single-cell protein data filtering for 120 dual time-point patients...")
    
    # Read sample ID mapping table
    mapping_file = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/process_bulk/sample_ids.tsv"
    print(f"Reading sample ID mapping: {mapping_file}")
    id_mapping = pd.read_csv(mapping_file, sep='\t')
    print(f"Samples in mapping table: {len(id_mapping)}")
    
    # Read core samples data (120 dual time-point patients)
    core_samples_file = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/sample_distribution/core_samples_10X_metabolomics_proteomics_20250704_173841.csv"
    print(f"Reading core samples data: {core_samples_file}")
    core_samples = pd.read_csv(core_samples_file)
    
    # Get target sample ID list (INCOV format)
    target_sample_ids = core_samples['Sample ID'].unique().tolist()
    
    # Convert sample ID format: INCOV{id}-1 -> INCOV{id}-BL, INCOV{id}-2 -> INCOV{id}-AC
    def convert_to_mapping_format(sample_id):
        if sample_id.endswith('-1'):
            return sample_id.replace('-1', '-BL')
        elif sample_id.endswith('-2'):
            return sample_id.replace('-2', '-AC')
        else:
            return sample_id
    
    converted_sample_ids = [convert_to_mapping_format(sid) for sid in target_sample_ids]
    
    # Find corresponding standardized sample_id (COVID_X_TX format) from mapping table
    target_mapping = id_mapping[id_mapping['original_filename'].isin(converted_sample_ids)]
    standardized_sample_ids = target_mapping['sample_id'].tolist()
    print(f"Mapped standardized sample IDs: {len(standardized_sample_ids)}")
    
    # Read complete single-cell protein data
    sc_pro_file = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/processed_omics_all/datablock_pro.tsv"
    print(f"\nReading single-cell protein data: {sc_pro_file}")
    sc_pro_data = pd.read_csv(sc_pro_file, sep='\t', index_col=0)
    print(f"Original single-cell protein data dimensions: {sc_pro_data.shape}")
    
    # Match samples using standardized IDs
    available_samples = sc_pro_data.columns.tolist()
    final_matched_samples = [s for s in standardized_sample_ids if s in available_samples]
    match_method = "standardized ID (COVID_X_TX)"
    
    print(f"\nSample matching results (using {match_method}):")
    print(f"  Target samples: {len(standardized_sample_ids)}")
    print(f"  Successfully matched: {len(final_matched_samples)}")
    
    # Filter data
    if final_matched_samples:
        filtered_data = sc_pro_data[final_matched_samples]
        print(f"\nFiltered data dimensions: {filtered_data.shape}")
        
        # Save filtered data to processed_omics_120 directory
        output_dir = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/processed_omics_120"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "sc_pro_120patients.csv")
        filtered_data.to_csv(output_file)
        print(f"\nFiltered single-cell protein data saved to: {output_file}")
        
        # Create sample metadata
        matched_core_sample_ids = []
        for std_id in final_matched_samples:
            orig_filename_row = target_mapping[target_mapping['sample_id'] == std_id]
            if not orig_filename_row.empty:
                orig_filename = orig_filename_row['original_filename'].iloc[0]
                if orig_filename.endswith('-BL'):
                    core_id = orig_filename.replace('-BL', '-1')
                elif orig_filename.endswith('-AC'):
                    core_id = orig_filename.replace('-AC', '-2')
                else:
                    core_id = orig_filename
                matched_core_sample_ids.append(core_id)
        
        sample_metadata = core_samples[core_samples['Sample ID'].isin(matched_core_sample_ids)].copy()
        
        metadata_file = os.path.join(output_dir, "sc_pro_120patients_metadata.csv")
        sample_metadata.to_csv(metadata_file, index=False)
        print(f"Sample metadata saved to: {metadata_file}")
        
        # Filter for patients with dual time points only
        patient_groups = sample_metadata.groupby('Study Subject ID').size()
        dual_timepoint_patients = patient_groups[patient_groups == 2].index.tolist()
        print(f"\nFound {len(dual_timepoint_patients)} dual time-point patients.")
        
        # Final filtering
        final_metadata = sample_metadata[sample_metadata['Study Subject ID'].isin(dual_timepoint_patients)]
        final_standardized_ids = []
        for sample_id in final_metadata['Sample ID']:
            converted_id = convert_to_mapping_format(sample_id)
            mapping_row = target_mapping[target_mapping['original_filename'] == converted_id]
            if not mapping_row.empty:
                final_standardized_ids.append(mapping_row['sample_id'].iloc[0])
        
        final_filtered_data = filtered_data[final_standardized_ids]
        print(f"Final filtered data dimensions: {final_filtered_data.shape}")

        # Save final files
        final_metadata.to_csv(metadata_file, index=False)
        final_filtered_data.to_csv(output_file)
        print(f"Final 120-patient data saved to files.")

        print(f"\n✅ Single-cell protein data filtering complete!")
        return output_file, metadata_file
    else:
        print("Error: No matching samples found!")
        return None, None

if __name__ == "__main__":
    main()