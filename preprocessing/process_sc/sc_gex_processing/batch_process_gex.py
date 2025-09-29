#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch processing for single-cell GEX data to save memory and disk space
This script coordinates the execution of process_gex_data.py
"""

import os
import glob
import subprocess
import shutil
from datetime import datetime

def log_message(message):
    """Print timestamped message"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def main():
    """Main scheduling function"""
    log_message("=== Starting GEX data batch processing ===")

    # Define paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir))))  # DIVAS-code directory
    WORKSPACE_ROOT = repo_root
    BASE_DIR = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing/process_sc/sc_gex_processing")
    
    SOURCE_DIR = os.path.join(BASE_DIR, "gex_data_unzipped")
    TEMP_DIR = os.path.join(BASE_DIR, "temp_gex_processing_batch")
    PROCESS_SCRIPT = os.path.join(BASE_DIR, "process_gex_data.py")
    
    # Check if source directory and script exist
    if not os.path.isdir(SOURCE_DIR):
        log_message(f"Error: Source directory '{SOURCE_DIR}' does not exist.")
        return
    if not os.path.isfile(PROCESS_SCRIPT):
        log_message(f"Error: Processing script '{PROCESS_SCRIPT}' does not exist.")
        return
        
    # Prepare batch processing
    all_files = sorted(glob.glob(os.path.join(SOURCE_DIR, "*.txt")))
    total_files = len(all_files)
    batch_size = 20
    num_batches = (total_files + batch_size - 1) // batch_size

    log_message(f"Found {total_files} files, will process in {num_batches} batches of {batch_size} each.")

    # Process each batch
    for i in range(num_batches):
        batch_num = i + 1
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, total_files)
        batch_files = all_files[start_idx:end_idx]

        log_message(f"Processing batch {batch_num}/{num_batches} ({len(batch_files)} files)")

        # Create temporary directory for this batch
        if os.path.exists(TEMP_DIR):
            shutil.rmtree(TEMP_DIR)
        os.makedirs(TEMP_DIR, exist_ok=True)

        try:
            # Copy batch files to temporary directory
            for file_path in batch_files:
                shutil.copy2(file_path, TEMP_DIR)
            
            log_message(f"Copied {len(batch_files)} files to temporary directory")

            # Run processing script on this batch
            cmd = ["python3", PROCESS_SCRIPT, TEMP_DIR]
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                log_message(f"Batch {batch_num} processed successfully")
            else:
                log_message(f"Error in batch {batch_num}: {result.stderr}")

        except Exception as e:
            log_message(f"Exception in batch {batch_num}: {str(e)}")

        finally:
            # Clean up temporary directory
            if os.path.exists(TEMP_DIR):
                shutil.rmtree(TEMP_DIR)
                log_message(f"Cleaned up temporary directory for batch {batch_num}")

    log_message("=== Batch processing complete ===")

if __name__ == "__main__":
    main()