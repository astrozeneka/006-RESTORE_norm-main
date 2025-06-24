import os

import scimap as sm
import pandas as pd
import anndata as ad
from glob import glob

# =============================================================================
# X. MAIN WORKFLOW FUNCTION
# =============================================================================
# =============================================================================
# MAIN WORKFLOW FUNCTION
# =============================================================================

def run_complete_phenotyping(data_path, phenotype_path="data/phenotype_workflow.csv"):
    """
    Run complete phenotyping workflow on IHC multiplex data

    Parameters:
    -----------
    data_path : str
        Path to the cell data CSV file
    phenotype_path : str
        Path to the phenotype workflow CSV file

    Returns:
    --------
    adata : AnnData
        Annotated data object with phenotyping results
    """

    try:
        # Load the data
        df = pd.read_csv(data_path)
        print(f"Loaded data with shape: {df.shape}")

        # Define feature columns (markers) - mapping to your new data structure
        feature_cols = [
            'CD3d',  # CD3
            'CD4',  # CD4
            'CD8',  # CD8
            'CD68',  # CD68
            'CD163',  # CD163
            'PD1',  # PD1
            'PanCK',  # PanCK
            'Ki67',  # Ki67
            'MHCII',  # MHC-II
            'CD11c',  # CD11c
            'CD20',  # CD20
            'CD31',  # CD31
            'FOXP3'  # FOXP3
        ]

        # Define metadata columns - mapping to your new data structure
        meta_cols = [
            'Object ID',  # CellID
            'Centroid X',  # X_centroid
            'Centroid Y',  # Y_centroid
            'Parent Region',  # Additional metadata
            'Parent Classification',
            'Classification',
            'Area'
        ]

        # Check which feature columns are actually present in the data
        available_features = [col for col in feature_cols if col in df.columns]
        missing_features = [col for col in feature_cols if col not in df.columns]

        print(f"Available feature columns: {len(available_features)}")
        print(f"Available features: {available_features}")
        if missing_features:
            print(f"Missing features: {missing_features}")

        # Check metadata columns
        available_meta = [col for col in meta_cols if col in df.columns]
        missing_meta = [col for col in meta_cols if col not in df.columns]

        print(f"Available metadata columns: {available_meta}")
        if missing_meta:
            print(f"Missing metadata: {missing_meta}")

        # Extract expression data and metadata
        expr = df[available_features].copy()
        meta = df[available_meta].copy()

        # Create simplified marker names for AnnData (remove ": Cell: Mean" suffix)
        marker_names = [col.replace(': Cell: Mean', '') for col in available_features]

        # Create AnnData object
        adata = ad.AnnData(expr.values)
        adata.var.index = marker_names  # simplified marker names

        # Set up observations (cells) metadata
        # Use Object ID as the index
        meta_indexed = meta.set_index('Object ID')
        adata.obs = meta_indexed

        # Rename coordinate columns to match scimap expectations
        if 'Centroid X' in adata.obs.columns:
            adata.obs['X_centroid'] = adata.obs['Centroid X']
        if 'Centroid Y' in adata.obs.columns:
            adata.obs['Y_centroid'] = adata.obs['Centroid Y']

        print(f"Created AnnData object with {adata.n_obs} cells and {adata.n_vars} markers")
        print(f"Markers: {list(adata.var.index)}")

        # Load phenotype workflow
        try:
            phenotype = pd.read_csv(phenotype_path)
            print(f"Loaded phenotype workflow with {len(phenotype)} phenotype definitions")

            # TODO: NORMALIZATION
            # TODO: add sm.pp.rescale(adata, method="quantile") # or method="standard" Here

            # You may need to update the phenotype workflow CSV to match your new marker names
            # The phenotype CSV should use the simplified marker names (without ": Cell: Mean")

            # Run phenotyping
            adata = sm.tl.phenotype_cells(
                adata=adata,
                phenotype=phenotype,
                # gate=0.5
            )

            # View results
            if 'phenotype' in adata.obs.columns:
                print("\nPhenotyping completed successfully!")
                print(f"Unique phenotypes found: {adata.obs['phenotype'].unique()}")
                print(f"Phenotype distribution:")
                print(adata.obs['phenotype'].value_counts())
            else:
                print("Warning: No 'phenotype' column found in results")

        except FileNotFoundError:
            print(f"Phenotype workflow file not found: {phenotype_path}")
            print("Returning AnnData object without phenotyping")
        except Exception as e:
            print(f"Error during phenotyping: {e}")
            print("Returning AnnData object without phenotyping")

        # extract basename from data_path
        import os
        base_name = os.path.basename(data_path).replace('.csv', '')
        slug = base_name.replace('cell_data_wih_regions_', '')
        # Svve to CSV
        adata.obs.to_csv(f'results/phenotyped_cells_{slug}.csv', index=True)

        return adata

    except Exception as e:
        print(f"Error in run_complete_phenotyping: {e}")
        return None


import argparse

# Load sample name from command
parser = argparse.ArgumentParser(description="Run complete phenotyping workflow")
parser.add_argument('--data-path', type=str, default=None,
                    help='Path to the cell data CSV file')
args = parser.parse_args()

def main(data_path):
    os.mkdir('results-phenotypes') if not os.path.exists('results-phenotypes') else None
    adata = run_complete_phenotyping(data_path)

    if adata is not None:
        print(f"\nFinal AnnData shape: {adata.shape}")
        print(f"Available observations (metadata): {list(adata.obs.columns)}")

        # Save results if needed
        slug = os.path.basename(data_path).replace('.csv', '')
        adata.write(f'results-phenotypes/phenotypes/{slug}.h5ad')

        # Example of accessing results
        if 'phenotype' in adata.obs.columns:
            print("\nSample of phenotyping results:")
            print(adata.obs[['X_centroid', 'Y_centroid', 'phenotype']].head(10))

        # Save data into file
        adata.obs.to_csv(f'results-phenotypes/phenotype_{slug}.csv', index=True)


    else:
        print("Failed to process data")

if __name__ == '__main__':
    os.mkdir("results-phenotypes/phenotypes")
    if args.data_path:
        main(args.data_path)
    else:
        data_path_list = glob("cell_data/normalized/*.csv")
        for data_path in data_path_list:
            print(f"Processing file: {data_path}")
            main(data_path)

