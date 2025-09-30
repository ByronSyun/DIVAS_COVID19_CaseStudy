import pandas as pd
import os

def verify_sample_order():
    """
    Verify sample order consistency across all omics data files
    """
    # Path configuration - script is now in the data directory
    script_dir = os.path.dirname(os.path.abspath(__file__))  # processed_omics_120 directory
    base_dir = script_dir
    
    files_to_check = {
        "GEX": {"path": os.path.join(base_dir, "sc_gex_120patients_aligned.csv"), "sep": ","},
        "Proteomics": {"path": os.path.join(base_dir, "proteomics_120patients.csv"), "sep": ","},
        "Metabolomics": {"path": os.path.join(base_dir, "metabolomics_120patients.csv"), "sep": ","},
        "sc-Proteomics": {"path": os.path.join(base_dir, "sc_pro_120patients.csv"), "sep": ","}
    }

    headers = {}
    
    print("Reading headers from all files...")
    all_files_found = True
    for name, info in files_to_check.items():
        try:
            # Read only header row to get column names
            df_header = pd.read_csv(info["path"], sep=info["sep"], engine='python', nrows=0)
            sample_ids = df_header.columns.tolist()[1:]
            headers[name] = sample_ids
            print(f"{name}: {len(sample_ids)} samples")
        except FileNotFoundError:
            print(f"Error: File not found for {name}: {info['path']}")
            all_files_found = False
        except Exception as e:
            print(f"Error reading {name}: {e}")
            all_files_found = False

    if not all_files_found or len(headers) < 2:
        print("\nVerification aborted: Not enough data to compare")
        return

    print("\nSample header summary:")
    for name, sample_ids in headers.items():
        print(f"{name}: {len(sample_ids)} samples")
        if len(sample_ids) > 10:
            print(f"  First 5: {sample_ids[:5]}")
            print(f"  Last 5: {sample_ids[-5:]}")
        
    print("\nComparing sample order...")
    
    # Use first header as reference
    reference_name = list(headers.keys())[0]
    reference_header = headers[reference_name]
    
    all_match = True
    for name_to_compare, header_to_compare in headers.items():
        if name_to_compare == reference_name:
            continue
        
        print(f"Comparing {reference_name} vs {name_to_compare}...")
        
        # Check length
        if len(reference_header) != len(header_to_compare):
            print(f"  Mismatch: Different sample counts")
            all_match = False
            continue

        # Check order
        if reference_header == header_to_compare:
            print(f"  Match: Sample order identical")
        else:
            print(f"  Mismatch: Sample order different")
            # Find first difference
            for i, (ref_id, comp_id) in enumerate(zip(reference_header, header_to_compare)):
                if ref_id != comp_id:
                    print(f"    First difference at index {i}: {ref_id} vs {comp_id}")
                    break
            all_match = False
            
    print("\nVerification result:")
    if all_match:
        print("Success: All files have identical sample order")
    else:
        print("Failure: Sample order is not consistent across files")
        print("Please use alignment script to fix order before DIVAS analysis")

if __name__ == "__main__":
    verify_sample_order() 