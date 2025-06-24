import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from glob import glob
import seaborn as sns


def extract_sample_id(file_path):
    """Extract sample ID from filename"""
    filename = Path(file_path).stem
    # Look for pattern like S19_12126B1 in the filename
    parts = filename.split('_')
    for i, part in enumerate(parts):
        if part.startswith('S') and len(part) > 1:
            # Try to find the sample ID pattern (e.g., S19_12126B1)
            if i + 1 < len(parts):
                sample_id = f"{part}_{parts[i + 1]}"
                return sample_id
            else:
                return part
    # Fallback: use the last part of filename before extension
    return parts[-1] if parts else "sample"


def plot_phenotype_spatial_distribution(data_path, output_dir="plots"):
    """
    Create and save a scatter plot of cell spatial distribution colored by phenotype
    """
    # Read the CSV file
    try:
        df = pd.read_csv(data_path)
        print(f"Loaded data with {len(df)} rows")
    except FileNotFoundError:
        print(f"Error: File {data_path} not found")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Check if required columns exist
    required_columns = ['X_centroid', 'Y_centroid', 'phenotype']
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        return

    # Remove rows with missing coordinates or phenotype
    initial_count = len(df)
    df = df.dropna(subset=['X_centroid', 'Y_centroid', 'phenotype'])
    final_count = len(df)

    if final_count < initial_count:
        print(f"Removed {initial_count - final_count} rows with missing data")

    if final_count == 0:
        print("Error: No valid data rows remaining after removing missing values")
        return

    # Get unique phenotypes and sort them alphabetically
    unique_phenotypes = sorted(df['phenotype'].unique())
    print(f"Found {len(unique_phenotypes)} unique phenotypes: {unique_phenotypes}")

    # Extract sample ID for title
    sample_id = extract_sample_id(data_path)

    # Create the scatter plot
    plt.figure(figsize=(12, 10))

    # Use recognizable colors for each phenotype
    recognizable_colors = [
        'red', 'blue', 'green', 'orange', 'purple', 'brown',
        'pink', 'gray', 'olive', 'cyan', 'magenta', 'yellow',
        'darkred', 'darkblue', 'darkgreen', 'darkorange', 'darkviolet',
        'maroon', 'navy', 'forestgreen', 'gold', 'indigo', 'crimson'
    ]

    # Plot each phenotype with a different color
    for i, phenotype in enumerate(unique_phenotypes):
        phenotype_data = df[df['phenotype'] == phenotype]
        color = recognizable_colors[i % len(recognizable_colors)]  # Cycle through colors if more phenotypes than colors

        plt.scatter(phenotype_data['X_centroid'],
                    phenotype_data['Y_centroid'],
                    c=color,
                    label=phenotype,
                    alpha=0.7,
                    s=2,  # Reduced point size
                    edgecolors='black',
                    linewidth=0.1)

        print(f"  {phenotype}: {len(phenotype_data)} cells")

    # Customize the plot
    plt.xlabel('X Coordinate', fontsize=12, fontweight='bold')
    plt.ylabel('Y Coordinate', fontsize=12, fontweight='bold')
    plt.title(f'Spatial Distribution of Cell Phenotypes - {sample_id}', fontsize=14, fontweight='bold')

    # Add legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, markerscale=10)

    # Add grid for better readability
    plt.grid(True, alpha=0.3, linestyle='--')

    # Set equal aspect ratio to preserve spatial relationships
    plt.axis('equal')

    # Adjust layout to prevent legend cutoff
    plt.tight_layout()

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Save the plot
    output_filename = f"phenotype_spatial_distribution_{sample_id}.png"
    output_path = os.path.join(output_dir, output_filename)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")

    # Optionally show the plot (comment out if running headless)
    # plt.show()

    plt.close()

    # Print summary statistics
    print(f"\nSpatial range:")
    print(f"  X: {df['X_centroid'].min():.2f} to {df['X_centroid'].max():.2f}")
    print(f"  Y: {df['Y_centroid'].min():.2f} to {df['Y_centroid'].max():.2f}")


def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Create spatial distribution plot of cell phenotypes")
    parser.add_argument('--data-path', type=str,
                        default=None,
                        help='Path to the cell data CSV file')
    parser.add_argument('--output-dir', type=str, default='results-phenotypes',
                        help='Directory to save the output plot (default: plots)')

    args = parser.parse_args()


    data_path = args.data_path
    if data_path:
        # Create the spatial distribution plot
        plot_phenotype_spatial_distribution(args.data_path, args.output_dir)
    else:
        data_path_list = glob("results-phenotypes/phenotype_*.csv")
        for data_path in data_path_list:
            plot_phenotype_spatial_distribution(data_path, args.output_dir)




if __name__ == "__main__":
    main()