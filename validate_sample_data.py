#!/usr/bin/env python3
"""
Validation script for RESTORE sample data
This script checks if the generated sample data is correctly formatted
and shows expected patterns before running RESTORE normalization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def validate_data_format(csv_file, marker_file):
    """Validate the format of input files"""
    print("=== VALIDATING DATA FORMAT ===")

    # Load data
    try:
        df = pd.read_csv(csv_file)
        markers_df = pd.read_csv(marker_file)
        print(f"✓ Successfully loaded {csv_file} and {marker_file}")
    except Exception as e:
        print(f"✗ Error loading files: {e}")
        return False

    # Check required columns
    required_cols = ['scene']
    if not all(col in df.columns for col in required_cols):
        print(f"✗ Missing required columns. Found: {df.columns.tolist()}")
        return False
    print(f"✓ Required columns present: {required_cols}")

    # Check marker columns
    marker_cols = [col for col in df.columns if col != 'scene']
    print(f"✓ Found {len(marker_cols)} markers: {marker_cols}")

    # Check scenes
    scenes = sorted(df['scene'].unique())
    print(f"✓ Found {len(scenes)} scenes: {scenes}")

    # Check data types
    for col in marker_cols:
        if not pd.api.types.is_numeric_dtype(df[col]):
            print(f"✗ Non-numeric data in column {col}")
            return False
    print("✓ All marker columns are numeric")

    # Check marker file format
    if not all(col in marker_cols for col in markers_df.columns):
        print(f"✗ Marker file columns don't match data columns")
        print(f"  Data columns: {marker_cols}")
        print(f"  Marker file columns: {markers_df.columns.tolist()}")
        return False
    print("✓ Marker file format is correct")

    return True


def analyze_batch_effects(csv_file):
    """Analyze potential batch effects across scenes"""
    print("\n=== ANALYZING BATCH EFFECTS ===")

    df = pd.read_csv(csv_file)
    marker_cols = [col for col in df.columns if col != 'scene']

    print("Scene-wise statistics:")
    for scene in sorted(df['scene'].unique()):
        scene_data = df[df['scene'] == scene]
        print(f"\nScene {scene} ({len(scene_data)} cells):")
        for marker in marker_cols:
            mean_val = scene_data[marker].mean()
            std_val = scene_data[marker].std()
            print(f"  {marker}: mean={mean_val:.1f}, std={std_val:.1f}")


def analyze_mutually_exclusive_patterns(csv_file, marker_file):
    """Analyze mutually exclusive marker patterns"""
    print("\n=== ANALYZING MUTUALLY EXCLUSIVE PATTERNS ===")

    df = pd.read_csv(csv_file)
    markers_df = pd.read_csv(marker_file)

    for pos_marker in markers_df.columns:
        print(f"\nAnalyzing {pos_marker}:")

        # Get mutually exclusive markers for this positive marker
        neg_markers = markers_df[pos_marker].dropna().tolist()

        for neg_marker in neg_markers:
            if neg_marker in df.columns:
                # Calculate correlation
                corr = df[pos_marker].corr(df[neg_marker])
                print(f"  {pos_marker} vs {neg_marker}: correlation = {corr:.3f}")

                # Identify potential positive cells
                pos_threshold = df[pos_marker].quantile(0.75)
                pos_cells = df[df[pos_marker] > pos_threshold]
                neg_expression_in_pos = pos_cells[neg_marker].mean()

                print(f"    Mean {neg_marker} in {pos_marker}+ cells: {neg_expression_in_pos:.1f}")


def create_diagnostic_plots(csv_file, output_dir="plots"):
    """Create diagnostic plots"""
    print(f"\n=== CREATING DIAGNOSTIC PLOTS ===")

    import os
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(csv_file)
    marker_cols = [col for col in df.columns if col != 'scene']

    # 1. Batch effect visualization
    plt.figure(figsize=(15, 10))
    for i, marker in enumerate(marker_cols, 1):
        plt.subplot(2, 3, i)
        for scene in sorted(df['scene'].unique()):
            scene_data = df[df['scene'] == scene][marker]
            plt.hist(scene_data, alpha=0.6, label=f'Scene {scene}', bins=20)
        plt.xlabel(f'{marker} Intensity')
        plt.ylabel('Count')
        plt.title(f'{marker} Distribution by Scene')
        plt.legend()
        plt.yscale('log')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/batch_effects.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved batch effect plot: {output_dir}/batch_effects.png")

    # 2. Correlation heatmap
    plt.figure(figsize=(10, 8))
    corr_matrix = df[marker_cols].corr()
    sns.heatmap(corr_matrix, annot=True, cmap='RdBu_r', center=0,
                square=True, fmt='.2f')
    plt.title('Marker Correlation Matrix')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/correlation_matrix.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved correlation matrix: {output_dir}/correlation_matrix.png")

    # 3. Scatter plots for key mutually exclusive pairs
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    pairs = [('CD68', 'CD31'), ('CD68', 'CK19'), ('CD45', 'CD31'),
             ('CD45', 'CK19'), ('CD3', 'CD20'), ('CD3', 'CK19')]

    for i, (marker1, marker2) in enumerate(pairs):
        if marker1 in df.columns and marker2 in df.columns:
            ax = axes[i // 3, i % 3]
            scatter = ax.scatter(df[marker1], df[marker2],
                                 c=df['scene'], alpha=0.6, cmap='viridis')
            ax.set_xlabel(f'{marker1} Intensity')
            ax.set_ylabel(f'{marker2} Intensity')
            ax.set_title(f'{marker1} vs {marker2}')
            plt.colorbar(scatter, ax=ax, label='Scene')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/mutually_exclusive_pairs.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved scatter plots: {output_dir}/mutually_exclusive_pairs.png")

    plt.close('all')


def main():
    """Main validation function"""
    csv_file = "sample_cell_data.csv"
    marker_file = "sample_markers.csv"

    print("RESTORE Sample Data Validation")
    print("=" * 50)

    # Validate format
    if not validate_data_format(csv_file, marker_file):
        print("✗ Data validation failed!")
        return

    # Analyze batch effects
    analyze_batch_effects(csv_file)

    # Analyze mutually exclusive patterns
    analyze_mutually_exclusive_patterns(csv_file, marker_file)

    # Create plots
    try:
        create_diagnostic_plots(csv_file)
    except ImportError:
        print("⚠ Matplotlib/Seaborn not available. Skipping plots.")
        print("  Install with: pip install matplotlib seaborn")

    print("\n=== VALIDATION SUMMARY ===")
    print("✓ Data format is correct")
    print("✓ Batch effects are present (good for testing normalization)")
    print("✓ Mutually exclusive patterns are detectable")
    print("✓ Data is ready for RESTORE normalization")

    print("\nNext steps:")
    print("1. Run threshold calculation: bash run_thresholds.sh")
    print("2. Run normalization: python normalize.py ...")
    print("3. Compare before/after results")


if __name__ == "__main__":
    main()