#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CellTypist annotation for scRNA-seq .h5ad files using the PBMC model.
Supports majority voting aggregation. Minimal logging; core logic unchanged.
"""

import os
import scanpy as sc
import celltypist
from celltypist import models
import sys
import argparse
import gc
import time
import numpy as np
import pandas as pd
import logging
import traceback

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('celltypist_annotation.log')
    ]
)
logger = logging.getLogger('celltypist_server')

try:
    import psutil
    def get_memory_usage():
        """Return the memory usage in MB."""
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss / 1024 / 1024
        return mem
except ImportError:
    logger.warning("psutil not installed. Memory usage tracking will be disabled.")
    def get_memory_usage():
        return 0

def annotate_data_with_celltypist(
    source_dir=None,
    output_dir=None,
    model_name='COVID19_Immune_Landscape.pkl',
    sample_to_process=None,
    use_majority_voting=True,
    save_memory=True
):
    """
    Annotate single-cell data using CellTypist with optional majority voting.
    
    Parameters:
        source_dir: Directory containing h5ad files
        output_dir: Directory for output files
        model_name: Name of CellTypist model to use
        sample_to_process: Optional; if provided, only process this specific sample
        use_majority_voting: Whether to use majority voting for cell type annotation
        save_memory: Whether to use memory-saving mode (recommended for large datasets)
    """
    start_time = time.time()
    
    # --- Configuration ---
    # The source directory for processed .h5ad files
    source_dir = '../../preprocessing/process_sc/sc_gex_processing/processed_gex_data'
    # The output directory for annotated files
    output_dir = 'annotated_data_majority_voting'
    # PBMC model
    model_name = 'Adult_COVID19_PBMC'

    # --- Setup ---
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Using majority_voting: {use_majority_voting}")

    # Load CellTypist model
    logger.info(f"Loading CellTypist model: {model_name}...")
    try:
        models.download_models(model=model_name)
        model = models.Model.load(model=model_name)
        logger.info("Model loaded")
    except Exception as e:
        logger.error(f"Error loading model: {e}")
        return

    # --- Processing Loop ---
    if sample_to_process:
        # Process only the specified sample
        sample_dirs = [sample_to_process]
        if not os.path.exists(os.path.join(source_dir, sample_to_process)):
            logger.error(f"Error: Sample directory {sample_to_process} not found in {source_dir}")
            return
    else:
        # Process all samples
        try:
            sample_dirs = [d for d in os.listdir(source_dir) 
                          if os.path.isdir(os.path.join(source_dir, d))]
        except Exception as e:
            logger.error(f"Error reading source directory: {e}")
            return
    
    logger.info(f"Found {len(sample_dirs)} sample directories to process")
    
    # Create summary data structure
    summary_data = []

    for sample in sample_dirs:
        sample_path = os.path.join(source_dir, sample)
        h5ad_file = f"{sample}_processed.h5ad"
        h5ad_path = os.path.join(sample_path, h5ad_file)

        if os.path.exists(h5ad_path):
            logger.info(f"\nProcessing sample {sample}...")
            try:
                # Load data
                logger.info(f"  Loading data: {h5ad_path}")
                adata = sc.read_h5ad(h5ad_path)
                logger.info(f"  Loaded {adata.n_obs} cells and {adata.n_vars} genes")

                # Predict cell types
                logger.info("  Running CellTypist prediction ...")
                predictions = celltypist.annotate(
                    adata, 
                    model=model, 
                    majority_voting=use_majority_voting
                )
                logger.info("  Prediction done")
                
                # Add prediction results to adata object
                adata.obs['celltypist_predicted_labels'] = predictions.predicted_labels['predicted_labels']
                
                # If using majority_voting, record these results
                if use_majority_voting:
                    adata.obs['celltypist_majority_voting'] = predictions.predicted_labels['majority_voting']
                    # Compatibility for downstream scripts expecting 'majority_voting'
                    adata.obs['majority_voting'] = predictions.predicted_labels['majority_voting']
                    # Record cell type distribution using majority voting
                    cell_counts = adata.obs['celltypist_majority_voting'].value_counts()
                    example_labels = adata.obs['celltypist_majority_voting'].unique()[:5].tolist()
                    
                    # Record majority_voting cell type distribution
                    for cell_type, count in cell_counts.items():
                        summary_data.append({
                            'sample': sample,
                            'cell_type': cell_type,
                            'count': count,
                            'percentage': count / adata.n_obs * 100
                        })
                else:
                    # Record cell type distribution using predicted labels
                    cell_counts = adata.obs['celltypist_predicted_labels'].value_counts()
                    example_labels = adata.obs['celltypist_predicted_labels'].unique()[:5].tolist()
                    
                    # Record predicted_labels cell type distribution
                    for cell_type, count in cell_counts.items():
                        summary_data.append({
                            'sample': sample,
                            'cell_type': cell_type,
                            'count': count,
                            'percentage': count / adata.n_obs * 100
                        })
                
                logger.info(f"  Annotation complete. Example labels: {example_labels}")

                # Save annotated data
                output_path = os.path.join(output_dir, f"{sample}_annotated.h5ad")
                logger.info(f"  Saving annotated data to: {output_path}")
                adata.write(output_path)
                logger.info("  Data saved")
                
                # Clean up memory
                del adata, predictions
                gc.collect()
                

            except Exception as e:
                logger.error(f"  Error processing {sample}: {e}")
                logger.error(traceback.format_exc())
        else:
            logger.warning(f"\nSkipping {sample}, no processed .h5ad file found")

    # Save summary table
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(output_dir, "celltypist_annotation_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        logger.info(f"Summary data saved to: {summary_path}")

    # Create cell type distribution visualization if using majority_voting
    if use_majority_voting and len(summary_data) > 0:
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            # Aggregate by cell type
            cell_type_summary = pd.DataFrame(summary_data).groupby('cell_type')['count'].sum().reset_index()
            cell_type_summary = cell_type_summary.sort_values('count', ascending=False)
            
            plt.figure(figsize=(12, 8))
            sns.barplot(x='count', y='cell_type', data=cell_type_summary.head(20))
            plt.title('Top 20 Cell Types (Majority Voting)')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'top_celltypes_majority_voting.png'), dpi=300)
            logger.info(f"Cell type distribution plot saved")
        except Exception as e:
            logger.warning(f"Visualization skipped: {e}")

    elapsed_time = time.time() - start_time
    logger.info(f"All samples processed in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate single-cell data with CellTypist (Server Version)')
    parser.add_argument('--source-dir', type=str, help='Source data directory path')
    parser.add_argument('--output-dir', type=str, help='Output directory path')
    parser.add_argument('--model', type=str, default='COVID19_Immune_Landscape.pkl', help='CellTypist model name')
    parser.add_argument('--sample', type=str, help='Process only this specific sample (e.g., INCOV018-BL)')
    parser.add_argument('--no-majority-voting', action='store_true', help='Disable majority voting')
    parser.add_argument('--save-memory', action='store_true', help='Enable memory-saving mode (recommended for large datasets)')
    args = parser.parse_args()
    
    annotate_data_with_celltypist(
        source_dir=args.source_dir,
        output_dir=args.output_dir,
        model_name=args.model,
        sample_to_process=args.sample,
        use_majority_voting=not args.no_majority_voting,
        save_memory=args.save_memory
    ) 