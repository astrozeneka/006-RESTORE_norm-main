from glob import glob
from os import path
import pandas as pd
import os


def merge_sampled_files():
    """
    Merge sampled CSV files into a single file for RESTORE_norm
    """
    # Get all sampled files
    sampled_files = glob("cell_data/sampled/*.csv")

    if not sampled_files:
        print("No sampled files found in csv_data/sampled/")
        return

    # Create output directory if it doesn't exist
    os.makedirs("cell_data", exist_ok=True)

    merged_data = []

    slug_ids = {}
    slug_curr_id = 1

    for file_path in sampled_files:
        # Extract slug from filename
        slug = path.basename(file_path).replace('.csv', '')
        if slug not in slug_ids:
            slug_ids[slug] = str(slug_curr_id)
            slug_curr_id += 1

        # Read the sampled file
        df = pd.read_csv(file_path)

        # Identify marker columns (exclude metadata and DAPI columns)
        metadata_cols = ['Object ID', 'Parent Region', 'Parent Classification',
                         'Centroid X', 'Centroid Y', 'Classification', 'Area']

        # Get all columns
        all_cols = df.columns.tolist()

        # Filter to get only marker columns
        marker_cols = []
        for col in all_cols:
            # Skip metadata columns
            if col in metadata_cols:
                continue
            # Skip DAPI columns
            if col.startswith('DAPI'):
                continue
            # Skip columns that are just dashes or empty
            if col in ['-', '--', '']:
                continue

            marker_cols.append(col)

        print(f"Processing {slug}: Found {len(marker_cols)} marker columns")
        print(f"  Marker columns: {marker_cols}")

        # Create scene column
        df_scene = df[marker_cols].copy()
        df_scene['scene'] = slug_ids[slug]

        # Reorder columns to have 'scene' first
        columns_ordered = ['scene'] + marker_cols
        df_scene = df_scene[columns_ordered]

        merged_data.append(df_scene)

        print(f"  Added {len(df_scene)} rows from {slug}")

    # Combine all dataframes
    final_df = pd.concat(merged_data, ignore_index=True)

    # Save merged file
    output_path = "cell_data/merged_samples.csv"
    final_df.to_csv(output_path, index=False)

    print(f"\nMerged file saved to: {output_path}")
    print(f"Total rows: {len(final_df)}")
    print(f"Total columns: {len(final_df.columns)}")
    print(f"Columns: {final_df.columns.tolist()}")
    print(f"Scenes: {final_df['scene'].unique()}")

    # Show sample of the data
    print(f"\nFirst 5 rows:")
    print(final_df.head())

    return final_df


if __name__ == '__main__':
    merge_sampled_files()