import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
from glob import glob

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


def plot_phenotype_histogram(data_path, output_dir="plots"):
    """
    Create and save a histogram of phenotype occurrences as percentages
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

    # Check if phenotype column exists
    if 'phenotype' not in df.columns:
        print("Error: 'phenotype' column not found in the CSV file")
        print(f"Available columns: {list(df.columns)}")
        return

    # Count phenotype occurrences
    phenotype_counts = df['phenotype'].value_counts()

    # Calculate percentages
    total_cells = len(df)
    phenotype_percentages = (phenotype_counts / total_cells * 100).round(2)

    # Sort alphabetically by phenotype name
    phenotype_percentages = phenotype_percentages.sort_index()

    print(f"Total cells: {total_cells}")
    print("Phenotype distribution:")
    for phenotype, percentage in phenotype_percentages.items():
        print(f"  {phenotype}: {percentage}% ({phenotype_counts[phenotype]} cells)")

    # Create the histogram
    plt.figure(figsize=(12, 8))
    bars = plt.bar(range(len(phenotype_percentages)),
                   phenotype_percentages.values,
                   alpha=0.7,
                   color='steelblue',
                   edgecolor='black',
                   linewidth=0.5)

    # Extract the sample Id
    sample_id = extract_sample_id(data_path)

    # Customize the plot
    plt.xlabel('Phenotype', fontsize=12, fontweight='bold')
    plt.ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    plt.title(f'Distribution of Cell Phenotypes - {sample_id}', fontsize=14, fontweight='bold')

    # Set x-axis labels
    plt.xticks(range(len(phenotype_percentages)),
               phenotype_percentages.index,
               rotation=45,
               ha='right')

    # Add percentage labels on top of bars
    for i, (phenotype, percentage) in enumerate(phenotype_percentages.items()):
        plt.text(i, percentage + 0.5, f'{percentage}%',
                 ha='center', va='bottom', fontweight='bold', fontsize=10)

    # Add grid for better readability
    plt.grid(axis='y', alpha=0.3, linestyle='--')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Extract sample ID and create output filename
    sample_id = extract_sample_id(data_path)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Save the plot
    output_filename = f"phenotype_histogram_{sample_id}.png"
    output_path = os.path.join(output_dir, output_filename)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")

    # Optionally show the plot (comment out if running headless)
    # plt.show()

    plt.close()


def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Create histogram of cell phenotype distribution")
    parser.add_argument('--data-path', type=str,
                        default=None,
                        help='Path to the cell data CSV file')
    parser.add_argument('--output-dir', type=str, default='results-phenotypes',
                        help='Directory to save the output plot (default: plots)')

    args = parser.parse_args()

    data_path =  args.data_path
    if data_path:
        # Create the histogram
        plot_phenotype_histogram(args.data_path, args.output_dir)
    else:
        data_path_list = glob('results-phenotypes/phenotype_*.csv')
        for data_path in data_path_list:
            plot_phenotype_histogram(data_path, args.output_dir)



if __name__ == "__main__":
    main()