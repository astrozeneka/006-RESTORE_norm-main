#!/usr/bin/env python3
"""
Split normalized data into separate files for each sample with cell metadata.
"""

import pandas as pd
from pathlib import Path

sample_ids = {
    "1": "S19_12126B1",
    "2": "S19_25142A1",
    "3": "S19_31776B1",
}


def split_normalized_data_with_metadata(normalized_file, plain_data_dir, output_dir):
    """
    Split normalized data by scene into separate CSV files with cell metadata.

    Args:
        normalized_file (str): Path to the input normalized data CSV
        plain_data_dir (str): Directory containing plain data files with metadata
        output_dir (str): Directory to save the split files
    """
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load the normalized data
    print(f"Loading normalized data from {normalized_file}...")
    normalized_data = pd.read_csv(normalized_file)
    print(f"Normalized data shape: {normalized_data.shape}")

    # Get marker columns (exclude scene column)
    marker_columns = [col for col in normalized_data.columns if col != 'scene']
    print(f"Markers: {marker_columns}")

    # Define metadata columns to extract from plain data
    metadata_columns = [
        'Object ID', 'Parent Region', 'Parent Classification',
        'Centroid X', 'Centroid Y', 'Classification', 'Area'
    ]

    # Process each sample
    for scene, sample_id in sample_ids.items():
        print(f"\nProcessing scene {scene} ({sample_id})...")

        # Filter normalized data for this scene
        scene_normalized = normalized_data[normalized_data['scene'] == int(scene)]

        if scene_normalized.empty:
            print(f"Warning: No normalized data found for scene {scene}")
            continue

        # Load corresponding plain data file
        plain_data_file = Path(plain_data_dir) / f"{sample_id}.csv"

        if not plain_data_file.exists():
            print(f"Warning: Plain data file not found: {plain_data_file}")
            continue

        print(f"Loading plain data from {plain_data_file}...")
        plain_data = pd.read_csv(plain_data_file)
        print(f"Plain data shape: {plain_data.shape}")

        # Verify row count matches
        if len(scene_normalized) != len(plain_data):
            print(f"Warning: Row count mismatch for scene {scene}")
            print(f"Normalized data: {len(scene_normalized)} rows")
            print(f"Plain data: {len(plain_data)} rows")
            continue

        # Extract metadata from plain data
        metadata = plain_data[metadata_columns].copy()

        # Extract normalized marker data (without scene column)
        markers_data = scene_normalized[marker_columns].reset_index(drop=True)

        # Combine metadata and normalized data
        combined_data = pd.concat([metadata, markers_data], axis=1)

        # Save to file
        output_file = Path(output_dir) / f"{sample_id}.csv"
        combined_data.to_csv(output_file, index=False)

        print(f"Saved {len(combined_data)} rows with {len(combined_data.columns)} columns to {output_file}")
        print(f"Columns: {list(combined_data.columns)}")


if __name__ == '__main__':
    # Input and output paths
    normalized_file = "results/merged_samples_RESTORE.csv"
    plain_data_dir = "cell_data/sampled"
    output_dir = "cell_data/normalized"

    try:
        split_normalized_data_with_metadata(normalized_file, plain_data_dir, output_dir)
        print("\nData splitting with metadata completed successfully!")

    except FileNotFoundError as e:
        print(f"Error: Could not find input file - {e}")
        print("Please ensure the normalized data file and plain data files exist at the specified paths.")
    except Exception as e:
        print(f"Error: {e}")