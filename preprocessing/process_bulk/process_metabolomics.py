#!/usr/bin/env python3
import os
import pandas as pd

# This script processes the metabolomics data and its associated metadata file.
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

def process_metabolomics_data():
    """
    Processes the main metabolomics data matrix.
    The column headers of this file are sample IDs that need standardization.
    """
    metabolomics_path = os.path.join(TARGET_DIR, "metabolomics_data.csv")
    if not os.path.exists(metabolomics_path):
        print(f"Error: Metabolomics data file not found at: {metabolomics_path}")
        return False

    metabolomics_df = pd.read_csv(metabolomics_path)

    # Standardize column headers (which are sample IDs)
    # The first column is assumed to be metabolite identifiers and is skipped.
    original_cols = metabolomics_df.columns.tolist()
    first_col_name = original_cols[0]
    sample_id_cols = original_cols[1:]
    
    new_sample_id_cols = [standardize_sample_id(col) for col in sample_id_cols]

    metabolomics_df.columns = [first_col_name] + new_sample_id_cols
    
    metabolomics_df.to_csv(metabolomics_path, index=False)
    print(f"Successfully processed and updated sample IDs in '{os.path.basename(metabolomics_path)}'.")
    return True

def process_metabolomics_metadata():
    """
    Processes the metabolomics metadata file.
    Standardizes the 'sample_id' column.
    """
    metadata_path = os.path.join(TARGET_DIR, "metabolomics_metadata.csv")
    if not os.path.exists(metadata_path):
        print(f"Warning: Metabolomics metadata file not found, skipping: {metadata_path}")
        return True # Not a fatal error

    metadata_df = pd.read_csv(metadata_path)
    
    if "sample_id" in metadata_df.columns:
        metadata_df["sample_id"] = metadata_df["sample_id"].apply(standardize_sample_id)
        metadata_df.to_csv(metadata_path, index=False)
        print("Successfully processed and updated sample IDs in 'metabolomics_metadata.csv'.")
    else:
        print("Warning: 'sample_id' column not found in 'metabolomics_metadata.csv'. Cannot update.")
    
    return True

if __name__ == "__main__":
    print("--- Processing Metabolomics Data and Metadata ---")
    
    data_success = process_metabolomics_data()
    metadata_success = process_metabolomics_metadata()
    
    if data_success and metadata_success:
        print("--- Metabolomics processing completed successfully. ---")
    else:
        print("--- Metabolomics processing failed. Please review the output. ---")
        exit(1) 