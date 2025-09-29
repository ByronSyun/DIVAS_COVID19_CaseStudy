import pandas as pd
import numpy as np
import os
from scipy.stats import gmean

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

def main():
    """
    Apply CLR normalization to single-cell proteomics data
    """
    # Path configuration for new git repo structure
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))  # DIVAS_COVID19_CaseStudy directory
    input_file = os.path.join(repo_root, "preprocessing", "processed_omics_120", "sc_pro_120patients.csv")
    output_file = os.path.join(repo_root, "preprocessing", "processed_omics_120", "sc_pro_120patients_clr_normalized.csv")

    print("Starting single-cell proteomics CLR normalization...")

    try:
        # Read input data
        df = pd.read_csv(input_file, index_col=0)
        print(f"Data loaded. Shape: {df.shape}")
        
        # Ensure index is named 'Feature'
        df.index.name = 'Feature'
        
        # Apply CLR transformation
        clr_df = clr_transform(df)
        
        # Save normalized data
        clr_df.to_csv(output_file)
        print(f"CLR-normalized data saved to: {output_file}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    main() 