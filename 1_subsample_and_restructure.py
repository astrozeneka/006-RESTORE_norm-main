from glob import glob
from os import path
import pandas as pd
import numpy as np
import os

cell_data_list = glob("cell_data/plain/*.csv")
# Slug begin with S19, S20, S21, etc. and end with .csv (not included)
# Extract basename
cell_data_slugs = [path.basename(f) for f in cell_data_list]

if __name__ == '__main__':
    # Create output directory if it doesn't exist
    os.makedirs("cell_data/sampled", exist_ok=True)

    for cell_data_path in cell_data_list:
        slug = path.basename(cell_data_path)
        slug = slug.replace('.csv', '')  # Remove the .csv extension
        slug = slug.replace('cell_data_with_regions_', '')

        # Read csv
        cell_data_df = pd.read_csv(cell_data_path)

        # Sample first 1000 rows (or all rows if less than 1000)
        N = 1000
        if N < 0:
            sampled_df = cell_data_df
        else:
            sampled_df = cell_data_df.head(N)

        # Rename columns - remove ": Cell: Mean" suffix
        new_columns = {}
        for col in sampled_df.columns:
            if ": Cell: Mean" in col:
                # Extract marker name by removing ": Cell: Mean"
                new_name = col.replace(": Cell: Mean", "")
                new_columns[col] = new_name
            # Keep other columns as they are

        # Apply column renaming
        sampled_df = sampled_df.rename(columns=new_columns)

        # Save to csv_data/sampled/
        output_path = f"cell_data/sampled/{slug}.csv"
        sampled_df.to_csv(output_path, index=False)

        print(f"Processed {slug}: {len(cell_data_df)} -> {len(sampled_df)} rows, saved to {output_path}")