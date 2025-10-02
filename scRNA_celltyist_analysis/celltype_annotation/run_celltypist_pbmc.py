#!/usr/bin/env python3

"""
CellTypist annotation for scRNA-seq .h5ad files using the PBMC model.
"""

import os
import scanpy as sc
import celltypist
from celltypist import models
import argparse
import gc
import numpy as np
import pandas as pd

def annotate_data_with_celltypist(
    source_dir=None,
    output_dir=None,
    model_name='COVID19_Immune_Landscape.pkl',
    sample_to_process=None,
    use_majority_voting=True,
    save_memory=True
):
    """Annotate single-cell data using CellTypist with majority voting."""
    
    # Configuration
    source_dir = '../../preprocessing/process_sc/sc_gex_processing/processed_gex_data'
    output_dir = 'annotated_data_majority_voting'
    model_name = 'Adult_COVID19_PBMC'

    # Setup
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Load CellTypist model
    print(f"Loading CellTypist model: {model_name}...")
    try:
        models.download_models(model=model_name)
        model = models.Model.load(model=model_name)
        print("Model loaded successfully")
    except Exception as e:
        print(f"Error loading model: {e}")
        return

    # --- Processing Loop ---
    if sample_to_process:
        # Process only the specified sample
        sample_dirs = [sample_to_process]
        if not os.path.exists(os.path.join(source_dir, sample_to_process)):
            print(f"Error: Sample directory {sample_to_process} not found in {source_dir}")
            return
    else:
        # Process all samples
        try:
            sample_dirs = [d for d in os.listdir(source_dir) 
                          if os.path.isdir(os.path.join(source_dir, d))]
        except Exception as e:
            print(f"Error reading source directory: {e}")
            return
    
    print(f"Found {len(sample_dirs)} sample directories to process")
    
    # Create summary data structure
    summary_data = []

    for sample in sample_dirs:
        sample_path = os.path.join(source_dir, sample)
        h5ad_file = f"{sample}_processed.h5ad"
        h5ad_path = os.path.join(sample_path, h5ad_file)

        if os.path.exists(h5ad_path):
            print(f"Processing sample {sample}...")
            try:
                # Load data
                adata = sc.read_h5ad(h5ad_path)
                print(f"  Loaded {adata.n_obs} cells and {adata.n_vars} genes")

                # Predict cell types
                predictions = celltypist.annotate(
                    adata, 
                    model=model, 
                    majority_voting=use_majority_voting
                )
                
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
                
                print(f"  Annotation complete. Example labels: {example_labels}")

                # Save annotated data
                output_path = os.path.join(output_dir, f"{sample}_annotated.h5ad")
                adata.write(output_path)
                print(f"  Saved to: {output_path}")
                
                # Clean up memory
                del adata, predictions
                gc.collect()

            except Exception as e:
                print(f"  Error processing {sample}: {e}")
        else:
            print(f"Skipping {sample}, no processed .h5ad file found")

    # Save summary table
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(output_dir, "celltypist_annotation_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"Summary data saved to: {summary_path}")

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
            print("Cell type distribution plot saved")
        except Exception as e:
            print(f"Visualization skipped: {e}")

    print("All samples processed successfully")

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