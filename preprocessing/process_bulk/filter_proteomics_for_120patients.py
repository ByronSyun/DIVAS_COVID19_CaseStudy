#!/usr/bin/env python3
"""
Filter proteomics data for 120 dual time-point patients
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

def convert_to_3digit_format(sample_id):
    """Convert standardized ID to 3-digit format used in proteomics data"""
    if 'COVID_' in sample_id:
        parts = sample_id.split('_')
        if len(parts) >= 3:
            try:
                patient_num = int(parts[1])
                time_point = parts[2]
                return f"COVID_{patient_num:03d}_{time_point}"
            except:
                return sample_id
    return sample_id

def main():
    print("Starting proteomics data filtering for 120 dual time-point patients...")
    
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
    print(f"Target samples: {len(target_sample_ids)}")
    print(f"Target sample ID examples: {target_sample_ids[:5]}")
    
    # Convert sample ID format: INCOV{id}-1 -> INCOV{id}-BL, INCOV{id}-2 -> INCOV{id}-AC
    def convert_to_mapping_format(sample_id):
        """Convert core sample format to mapping table format"""
        if sample_id.endswith('-1'):
            return sample_id.replace('-1', '-BL')
        elif sample_id.endswith('-2'):
            return sample_id.replace('-2', '-AC')
        else:
            return sample_id  # Keep original format
    
    # Convert to mapping format
    converted_sample_ids = [convert_to_mapping_format(sid) for sid in target_sample_ids]
    print(f"Converted sample ID examples: {converted_sample_ids[:5]}")
    
    # Find corresponding standardized sample_id (COVID_X_TX format) from mapping table
    target_mapping = id_mapping[id_mapping['original_filename'].isin(converted_sample_ids)]
    standardized_sample_ids = target_mapping['sample_id'].tolist()
    print(f"Mapped standardized sample IDs: {len(standardized_sample_ids)}")
    print(f"Standardized sample ID examples: {standardized_sample_ids[:5]}")
    
    # Convert to 3-digit format for proteomics data
    proteomics_format_ids = [convert_to_3digit_format(sid) for sid in standardized_sample_ids]
    print(f"Proteomics format sample ID examples: {proteomics_format_ids[:5]}")
    
    # Read complete proteomics data
    proteomics_file = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/processed_omics_all/improved_proteomics_data.csv"
    print(f"\nReading proteomics data: {proteomics_file}")
    proteomics_data = pd.read_csv(proteomics_file, index_col=0)
    print(f"Original proteomics data dimensions: {proteomics_data.shape}")
    print(f"Proteomics sample ID examples: {proteomics_data.columns[:5].tolist()}")
    
    # Check sample ID formats in proteomics data
    available_samples = proteomics_data.columns.tolist()
    
    print("\n=== Sample ID matching analysis ===")
    
    # Method 1: Use converted 3-digit format matching
    matched_samples = [s for s in proteomics_format_ids if s in available_samples]
    print(f"3-digit format matching: {len(matched_samples)}")
    
    # Method 2: Direct matching with standardized ID (COVID_X_TX)
    matched_standard = [s for s in standardized_sample_ids if s in available_samples]
    print(f"Standardized ID direct matching: {len(matched_standard)}")
    
    # Choose best matching method
    if len(matched_samples) >= len(matched_standard):
        final_matched_samples = matched_samples
        match_method = "3-digit format (COVID_001_TX)"
        # Create mapping: proteomics format -> standardized format
        id_conversion_map = dict(zip(proteomics_format_ids, standardized_sample_ids))
    else:
        final_matched_samples = matched_standard
        match_method = "standardized ID (COVID_X_TX)"
        id_conversion_map = dict(zip(standardized_sample_ids, standardized_sample_ids))
    
    missing_samples = [s for s in proteomics_format_ids if s not in available_samples]
    
    print(f"\nSample matching results (using {match_method}):")
    print(f"  Target samples: {len(proteomics_format_ids)}")
    print(f"  Successfully matched: {len(final_matched_samples)}")
    print(f"  Missing samples: {len(missing_samples)}")
    
    if missing_samples and len(missing_samples) <= 20:
        print(f"Missing samples: {missing_samples}")
    elif missing_samples:
        print(f"Missing sample examples (first 10): {missing_samples[:10]}")
    
    # Filter proteomics data
    if final_matched_samples:
        filtered_proteomics = proteomics_data[final_matched_samples]
        print(f"\nFiltered proteomics data dimensions: {filtered_proteomics.shape}")
        
        # Convert column names back to standardized format for consistency
        if match_method.startswith("3-digit format"):
            # Rename columns to standardized format
            column_rename_map = {}
            for prot_id in final_matched_samples:
                if prot_id in id_conversion_map:
                    column_rename_map[prot_id] = id_conversion_map[prot_id]
            
            filtered_proteomics = filtered_proteomics.rename(columns=column_rename_map)
            print(f"Column names converted to standardized format")
        
        # Check data quality
        missing_values_per_sample = filtered_proteomics.isnull().sum()
        missing_values_per_protein = filtered_proteomics.isnull().sum(axis=1)
        
        print(f"\nData quality check:")
        print(f"  Sample missing values: min={missing_values_per_sample.min()}, max={missing_values_per_sample.max()}, mean={missing_values_per_sample.mean():.2f}")
        print(f"  Protein missing values: min={missing_values_per_protein.min()}, max={missing_values_per_protein.max()}, mean={missing_values_per_protein.mean():.2f}")
        
        # Save filtered data to processed_omics_120 directory
        output_dir = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/processed_omics_120"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "proteomics_120patients.csv")
        filtered_proteomics.to_csv(output_file)
        print(f"\nFiltered proteomics data saved to: {output_file}")
        
        # Create sample metadata
        # Need to map matched standardized IDs back to core sample format
        matched_standardized_ids = list(filtered_proteomics.columns)
        matched_core_sample_ids = []
        
        for std_id in matched_standardized_ids:
            # Find corresponding original_filename from target_mapping
            orig_filename_row = target_mapping[target_mapping['sample_id'] == std_id]
            if len(orig_filename_row) > 0:
                orig_filename = orig_filename_row['original_filename'].iloc[0]
                # Convert back to core sample format
                if orig_filename.endswith('-BL'):
                    core_id = orig_filename.replace('-BL', '-1')
                elif orig_filename.endswith('-AC'):
                    core_id = orig_filename.replace('-AC', '-2')
                else:
                    core_id = orig_filename
                matched_core_sample_ids.append(core_id)
        
        sample_metadata = core_samples[core_samples['Sample ID'].isin(matched_core_sample_ids)].copy()
        
        metadata_file = os.path.join(output_dir, "proteomics_120patients_metadata.csv")
        sample_metadata.to_csv(metadata_file, index=False)
        print(f"Sample metadata saved to: {metadata_file}")
        
        # Generate summary report
        print(f"\n=== Proteomics data filtering summary ===")
        print(f"Matching method: {match_method}")
        print(f"Original proteomics data: {proteomics_data.shape[0]} proteins × {proteomics_data.shape[1]} samples")
        print(f"Filtered data: {filtered_proteomics.shape[0]} proteins × {filtered_proteomics.shape[1]} samples")
        print(f"Data retention rate: {len(final_matched_samples)}/{len(proteomics_format_ids)} ({len(final_matched_samples)/len(proteomics_format_ids)*100:.1f}%)")
        
        # Group by patient statistics
        if len(sample_metadata) > 0:
            patient_groups = sample_metadata.groupby('Study Subject ID').size()
            
            # Filter for patients with dual time points only
            dual_timepoint_patients = patient_groups[patient_groups == 2].index.tolist()
            print(f"\nFiltering for patients with dual time points only...")
            print(f"  Found {len(dual_timepoint_patients)} dual time-point patients.")
            
            # Filter metadata and proteomics data
            final_metadata = sample_metadata[sample_metadata['Study Subject ID'].isin(dual_timepoint_patients)]
            
            # Get standardized sample IDs from final metadata
            final_standardized_ids = []
            for sample_id in final_metadata['Sample ID']:
                converted_id = convert_to_mapping_format(sample_id)
                mapping_row = target_mapping[target_mapping['original_filename'] == converted_id]
                if not mapping_row.empty:
                    final_standardized_ids.append(mapping_row['sample_id'].iloc[0])

            # Select these samples from original filtered proteomics data
            final_filtered_proteomics = filtered_proteomics[final_standardized_ids]

            print(f"Final filtered data dimensions: {final_filtered_proteomics.shape}")

            # Re-save files
            final_metadata.to_csv(metadata_file, index=False)
            final_filtered_proteomics.to_csv(output_file)
            print(f"Final 120-patient data re-saved to files.")

            patient_groups = final_metadata.groupby('Study Subject ID').size()
            print(f"\nFinal patient time point distribution:")
            print(f"  Patients with 2 time points: {(patient_groups == 2).sum()}")
            print(f"  Patients with 1 time point: {(patient_groups == 1).sum()}")
            print(f"  Total patients: {len(patient_groups)}")
            
            # Statistics by time point
            timepoint_counts = final_metadata['Blood draw time point'].value_counts()
            print(f"\nTime point distribution:")
            for tp, count in timepoint_counts.items():
                print(f"  {tp}: {count} samples")
        
        return output_file, metadata_file
    else:
        print("Error: No matching samples found!")
        return None, None

if __name__ == "__main__":
    output_file, metadata_file = main()
    if output_file:
        print(f"\n✅ Proteomics data filtering complete!")
        print(f"📁 Data file: {output_file}")
        print(f"📁 Metadata file: {metadata_file}")
    else:
        print("❌ Proteomics data filtering failed!")