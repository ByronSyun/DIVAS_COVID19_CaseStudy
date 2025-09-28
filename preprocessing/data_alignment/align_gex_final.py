import pandas as pd
import os

def align_gex_to_standard():
    """
    Align GEX data to match sample order of other omics files
    """
    # Path configuration for new git repo structure
    base_dir = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/processed_omics_120"
    
    # Reference file with correct sample order
    reference_file = os.path.join(base_dir, "proteomics_120patients.csv")
    
    # GEX file to be aligned (from processed_omics_all)
    gex_input_dir = "/Users/byronsun/Desktop/DIVAS-code/DIVAS_COVID19_CaseStudy/preprocessing/processed_omics_all"
    gex_file_to_fix = os.path.join(gex_input_dir, "gex_divas_datablock.tsv")
    
    # Output aligned GEX file
    output_gex_file = os.path.join(base_dir, "sc_gex_120patients_aligned.csv")

    print("Starting GEX data alignment...")

    try:
        # Get standard sample order from reference file
        df_ref_header = pd.read_csv(reference_file, nrows=0)
        standard_order = df_ref_header.columns.tolist()[1:]
        print(f"Standard sample order: {len(standard_order)} samples")

        # Read GEX file
        gex_df = pd.read_csv(gex_file_to_fix, sep='\\t', index_col=0, engine='python')
        print(f"GEX data loaded. Shape: {gex_df.shape}")

        # Align GEX data to standard order
        aligned_gex_df = gex_df[standard_order]
        print(f"GEX data aligned. New shape: {aligned_gex_df.shape}")
        
        # Set feature column name
        aligned_gex_df.index.name = 'Feature'

        # Save aligned GEX file
        aligned_gex_df.to_csv(output_gex_file, sep=',')
        print(f"Aligned GEX data saved to: {output_gex_file}")

    except FileNotFoundError as e:
        print(f"Error: Required file not found. {e}")
    except KeyError as e:
        print(f"Error: Sample ID mismatch. {e}")
    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    align_gex_to_standard() 