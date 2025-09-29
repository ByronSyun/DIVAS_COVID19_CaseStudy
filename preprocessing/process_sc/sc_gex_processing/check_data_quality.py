#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Check GEX datablock data quality
Analyzes missing values and zero value proportions, provides data quality statistics
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Path configuration
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir))))  # DIVAS-code directory
WORKSPACE_ROOT = repo_root
BASE_DIR = os.path.join(WORKSPACE_ROOT, "DIVAS_COVID19_CaseStudy/preprocessing/process_sc/sc_gex_processing")
DATABLOCK_FILE = os.path.join(BASE_DIR, "gex_divas_datablock.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "quality_reports")

os.makedirs(OUTPUT_DIR, exist_ok=True)
log_file = os.path.join(OUTPUT_DIR, "data_quality_check.log")

def log_message(message):
    """Log message to file and console"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)

def check_datablock_quality(datablock_file):
    """Check data quality of datablock file"""
    log_message(f"Reading datablock file: {datablock_file}")
    
    # Read data
    try:
        data = pd.read_csv(datablock_file, sep='\t', index_col=0)
        log_message(f"Successfully read data, dimensions: {data.shape} (genes x samples)")
    except Exception as e:
        log_message(f"Error reading data: {str(e)}")
        return
    
    # Check missing values
    null_count = data.isnull().sum().sum()
    null_percentage = null_count / data.size * 100
    log_message(f"Missing values: {null_count} ({null_percentage:.2f}%)")
    
    # Check zero values
    zero_count = (data == 0).sum().sum()
    zero_percentage = zero_count / data.size * 100
    log_message(f"Zero values: {zero_count} ({zero_percentage:.2f}%)")
    
    # Basic statistics
    log_message(f"Total genes: {data.shape[0]}")
    log_message(f"Total samples: {data.shape[1]}")
    
    # Gene statistics
    gene_stats = {}
    gene_stats['genes_with_all_zeros'] = (data == 0).all(axis=1).sum()
    gene_stats['genes_with_some_expression'] = (data > 0).any(axis=1).sum()
    gene_stats['avg_expression_per_gene'] = data.mean(axis=1).mean()
    gene_stats['median_expression_per_gene'] = data.mean(axis=1).median()
    
    log_message(f"Genes with all zeros: {gene_stats['genes_with_all_zeros']}")
    log_message(f"Genes with some expression: {gene_stats['genes_with_some_expression']}")
    log_message(f"Average expression per gene: {gene_stats['avg_expression_per_gene']:.4f}")
    log_message(f"Median expression per gene: {gene_stats['median_expression_per_gene']:.4f}")
    
    # Sample statistics
    sample_stats = {}
    sample_stats['avg_expression_per_sample'] = data.mean(axis=0).mean()
    sample_stats['median_expression_per_sample'] = data.mean(axis=0).median()
    sample_stats['expressed_genes_per_sample'] = (data > 0).sum(axis=0).mean()
    
    log_message(f"Average expression per sample: {sample_stats['avg_expression_per_sample']:.4f}")
    log_message(f"Median expression per sample: {sample_stats['median_expression_per_sample']:.4f}")
    log_message(f"Average expressed genes per sample: {sample_stats['expressed_genes_per_sample']:.1f}")
    
    # Create visualizations
    create_quality_plots(data)
    
    # Save quality report
    save_quality_report(data, gene_stats, sample_stats)
    
    log_message("Data quality check complete")

def create_quality_plots(data):
    """Create quality control plots"""
    log_message("Creating quality control plots...")
    
    # Plot 1: Distribution of expression values
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    non_zero_values = data.values[data.values > 0]
    plt.hist(np.log10(non_zero_values + 1), bins=50, alpha=0.7)
    plt.xlabel('Log10(Expression + 1)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Non-zero Expression Values')
    
    # Plot 2: Number of expressed genes per sample
    plt.subplot(1, 3, 2)
    expressed_genes = (data > 0).sum(axis=0)
    plt.hist(expressed_genes, bins=30, alpha=0.7)
    plt.xlabel('Number of Expressed Genes')
    plt.ylabel('Number of Samples')
    plt.title('Expressed Genes per Sample')
    
    # Plot 3: Expression level per sample
    plt.subplot(1, 3, 3)
    sample_means = data.mean(axis=0)
    plt.hist(sample_means, bins=30, alpha=0.7)
    plt.xlabel('Mean Expression Level')
    plt.ylabel('Number of Samples')
    plt.title('Mean Expression per Sample')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "expression_distributions.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 4: Heatmap of top variable genes
    gene_var = data.var(axis=1)
    top_var_genes = gene_var.nlargest(50).index
    
    plt.figure(figsize=(15, 8))
    sns.heatmap(data.loc[top_var_genes], cmap='viridis', robust=True)
    plt.title('Heatmap: Top 50 Most Variable Genes')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.savefig(os.path.join(OUTPUT_DIR, "top_variable_genes_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    log_message("Quality control plots saved")

def save_quality_report(data, gene_stats, sample_stats):
    """Save comprehensive quality report"""
    log_message("Saving quality report...")
    
    report = {
        'Data_Dimensions': {
            'Genes': data.shape[0],
            'Samples': data.shape[1],
            'Total_Values': data.size
        },
        'Missing_Values': {
            'Count': data.isnull().sum().sum(),
            'Percentage': data.isnull().sum().sum() / data.size * 100
        },
        'Zero_Values': {
            'Count': (data == 0).sum().sum(),
            'Percentage': (data == 0).sum().sum() / data.size * 100
        },
        'Gene_Statistics': gene_stats,
        'Sample_Statistics': sample_stats
    }
    
    # Save as JSON for easy reading
    import json
    report_file = os.path.join(OUTPUT_DIR, "quality_report.json")
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Also save as text for human reading
    text_report_file = os.path.join(OUTPUT_DIR, "quality_report.txt")
    with open(text_report_file, 'w') as f:
        f.write("GEX Data Quality Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Data Dimensions:\n")
        f.write(f"  Genes: {data.shape[0]}\n")
        f.write(f"  Samples: {data.shape[1]}\n")
        f.write(f"  Total values: {data.size}\n\n")
        
        f.write(f"Missing Values:\n")
        f.write(f"  Count: {data.isnull().sum().sum()}\n")
        f.write(f"  Percentage: {data.isnull().sum().sum() / data.size * 100:.2f}%\n\n")
        
        f.write(f"Zero Values:\n")
        f.write(f"  Count: {(data == 0).sum().sum()}\n")
        f.write(f"  Percentage: {(data == 0).sum().sum() / data.size * 100:.2f}%\n\n")
        
        f.write(f"Gene Statistics:\n")
        for key, value in gene_stats.items():
            f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        f.write(f"Sample Statistics:\n")
        for key, value in sample_stats.items():
            f.write(f"  {key}: {value}\n")
    
    log_message(f"Quality report saved to: {report_file}")
    log_message(f"Text report saved to: {text_report_file}")

def main():
    """Main function"""
    log_message("Starting GEX data quality check")
    
    if not os.path.exists(DATABLOCK_FILE):
        log_message(f"Error: Datablock file not found: {DATABLOCK_FILE}")
        log_message("Please run create_divas_input.py first")
        return
    
    check_datablock_quality(DATABLOCK_FILE)
    
    log_message("Data quality check completed")

if __name__ == "__main__":
    main()