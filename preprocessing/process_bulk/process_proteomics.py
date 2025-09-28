#!/usr/bin/env python3
import os
import pandas as pd

# This script processes the proteomics data and associated metadata files.
# Its main function is to standardize sample IDs from the source format
# (e.g., 'INCOV001-BL') to a consistent format ('COVID_1_T1').

# Get the directory where the script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TARGET_DIR = SCRIPT_DIR

def standardize_sample_id(sample_id):
    """
    Standardizes a single sample ID.
    Converts 'INCOV...' format to 'COVID_...' format.
    Example: INCOV001-BL -> COVID_1_T1
    """
    if isinstance(sample_id, str) and sample_id.startswith("INCOV"):
        parts = sample_id.split("-")
        if len(parts) == 2:
            # Extract number after 'INCOV' and remove leading zeros
            sample_num = str(int(parts[0][5:]))
            time_point_code = parts[1]
            
            # Convert timepoint codes
            if time_point_code == "BL":
                time_point = "T1"
            elif time_point_code == "AC":
                time_point = "T2"
            else:
                time_point = time_point_code # Fallback
            return f"COVID_{sample_num}_{time_point}"
    return sample_id

def process_proteomics_data():
    """
    Processes the main proteomics data matrix.
    The column headers of this file are sample IDs that need standardization.
    """
    proteomics_path = os.path.join(TARGET_DIR, "proteomics_data.csv")
    if not os.path.exists(proteomics_path):
        print(f"Error: Proteomics data file not found at: {proteomics_path}")
        return False

    proteomics_df = pd.read_csv(proteomics_path)
    
    # Standardize column headers (which are sample IDs)
    # The first column is assumed to be protein identifiers and is skipped.
    original_cols = proteomics_df.columns.tolist()
    first_col_name = original_cols[0]
    sample_id_cols = original_cols[1:]
    
    new_sample_id_cols = [standardize_sample_id(col) for col in sample_id_cols]
    
    proteomics_df.columns = [first_col_name] + new_sample_id_cols
    
    proteomics_df.to_csv(proteomics_path, index=False)
    print(f"Successfully processed and updated sample IDs in '{os.path.basename(proteomics_path)}'.")
    return True

def process_metadata():
    """
    Processes all metadata files.
    Standardizes the 'sample_id' column in each file.
    """
    metadata_files = [
        "proteomics_metadata.csv",
        "clinical_metadata.csv",
        "healthy_donor_metadata.csv"
    ]
    
    all_successful = True
    for filename in metadata_files:
        metadata_path = os.path.join(TARGET_DIR, filename)
        if not os.path.exists(metadata_path):
            print(f"Warning: Metadata file not found, skipping: {metadata_path}")
            continue
            
        metadata_df = pd.read_csv(metadata_path)
        
        if "sample_id" in metadata_df.columns:
            metadata_df["sample_id"] = metadata_df["sample_id"].apply(standardize_sample_id)
            metadata_df.to_csv(metadata_path, index=False)
            print(f"Successfully processed and updated sample IDs in '{filename}'.")
        else:
            print(f"Warning: 'sample_id' column not found in '{filename}'. Cannot update.")
            all_successful = False

    return all_successful

if __name__ == "__main__":
    print("--- Processing Proteomics and Metadata Files ---")
    
    proteomics_success = process_proteomics_data()
    metadata_success = process_metadata()
    
    if proteomics_success and metadata_success:
        print("--- Proteomics and metadata processing completed successfully. ---")
    else:
        print("--- Proteomics and metadata processing completed with warnings. Please review the output. ---")
        exit(1) 