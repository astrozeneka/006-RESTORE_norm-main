#!/usr/bin/env python3
"""
Intensity Distribution Plotter for Cytometry Data
Generates distribution plots for raw and normalized data across samples and markers.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import argparse


def load_data(filepath):
    """Load CSV data and return DataFrame."""
    return pd.read_csv(filepath)


def create_distribution_plot(data, title_suffix, output_path):
    """
    Create intensity distribution plots for all markers across samples.

    Args:
        data (pd.DataFrame): Input data with scene and marker columns
        title_suffix (str): Suffix for plot title (e.g., "Raw" or "Normalized")
        output_path (str): Path to save the output plot
    """
    # Get marker columns (exclude scene column)
    markers = [col for col in data.columns if col != 'scene']

    # Get unique samples
    samples = sorted(data['scene'].unique())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green

    # Calculate subplot grid dimensions
    n_markers = len(markers)
    n_cols = 4
    n_rows = (n_markers + n_cols - 1) // n_cols

    # Create figure and subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4 * n_rows))
    fig.suptitle(f'Intensity Distributions - {title_suffix} Data', fontsize=16, y=0.98)

    # Flatten axes array for easier indexing
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes.flatten()
    else:
        axes = axes.flatten()

    # Plot distribution for each marker
    for i, marker in enumerate(markers):
        ax = axes[i]

        # Plot distribution for each sample
        for j, sample in enumerate(samples):
            sample_data = data[data['scene'] == sample][marker]

            # Plot KDE
            sample_data.plot.density(
                ax=ax,
                color=colors[j],
                alpha=0.7,
                linewidth=2,
                label=f'Sample {sample}'
            )

        ax.set_title(marker, fontsize=12, fontweight='bold')
        ax.set_xlabel('Intensity')
        ax.set_ylabel('Density')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for i in range(n_markers, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to: {output_path}")


def main():
    """Main function to generate both raw and normalized distribution plots."""

    # Create output directory
    Path("plots").mkdir(exist_ok=True)

    # File paths
    raw_data_path = "cell_data/merged_samples.csv"
    normalized_data_path = "results/merged_samples_RESTORE.csv"  # Assuming this is the normalized data path

    try:
        # Load raw data
        print("Loading raw data...")
        raw_data = load_data(raw_data_path)
        print(f"Raw data shape: {raw_data.shape}")

        # Create raw data distribution plot
        print("Generating raw data distribution plot...")
        create_distribution_plot(
            raw_data,
            "Raw",
            "plots/intensity_distribution_raw.png"
        )

        # Load normalized data
        print("Loading normalized data...")
        normalized_data = load_data(normalized_data_path)
        print(f"Normalized data shape: {normalized_data.shape}")

        # Create normalized data distribution plot
        print("Generating normalized data distribution plot...")
        create_distribution_plot(
            normalized_data,
            "Normalized",
            "plots/intensity_distribution_normalized.png"
        )

        print("All plots generated successfully!")

    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}")
        print("Please ensure the data files exist in the specified paths.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()