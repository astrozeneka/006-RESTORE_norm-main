#!/usr/bin/env python3
"""
Complete SCIMAP Analysis Workflow for IHC Data - ADAPTED VERSION
===============================================================
Author: Your Name
Date: June 2025

This script is adapted for your specific CSV format with QuPath-exported data.
Includes quality control, phenotyping, spatial analysis, and statistics.
"""

# =============================================================================
# 1. IMPORT LIBRARIES AND SETUP
# =============================================================================

import pandas as pd
import numpy as np
import scanpy as sc
import scimap as sm
import matplotlib.pyplot as plt
import seaborn as sns
from  glob import glob
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set up scanpy settings
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Set working directory
import os
# os.chdir('/path/to/your/data')  # Uncomment and set your path

print("✅ Libraries loaded successfully!")


# =============================================================================
# 2. DATA LOADING AND PREPROCESSING - ADAPTED FOR YOUR FORMAT
# =============================================================================

def load_and_preprocess_data(file_path):
    """
    Load and preprocess your QuPath-exported CSV data

    Parameters:
    -----------
    file_path : str
        Path to your CSV file containing single-cell data
    """

    print("📁 Loading and preprocessing data...")
    df = pd.read_csv(file_path)

    # Display original columns for reference
    print("📋 Original columns in your data:")
    for i, col in enumerate(df.columns):
        print(f"   {i + 1:2d}. {col}")

    # Step 1: Rename coordinate columns to scimap standard
    df = df.rename(columns={
        'Centroid X': 'X_centroid',
        'Centroid Y': 'Y_centroid',
        'Parent Region': 'imageid',
        'Object ID': 'cell_id'
    })

    # Step 2: Clean up marker names and select only the relevant ones
    marker_mapping = {
        'Ki67: Cell: Mean': 'Ki67',
        'CD11c: Cell: Mean': 'CD11c',
        'CD3d: Cell: Mean': 'CD3',
        'MHCII: Cell: Mean': 'MHC-II',
        'CD68: Cell: Mean': 'CD68',
        'CD8: Cell: Mean': 'CD8',
        'PD1: Cell: Mean': 'PD1',
        'FOXP3: Cell: Mean': 'FOXP3',
        'CD4: Cell: Mean': 'CD4',
        'CD20: Cell: Mean': 'CD20',
        'PanCK: Cell: Mean': 'PanCK',
        'CD163: Cell: Mean': 'CD163',
        'CD31: Cell: Mean': 'CD31'
    }

    # Rename marker columns
    df = df.rename(columns=marker_mapping)

    # Step 3: Select relevant columns for analysis
    coordinate_cols = ['X_centroid', 'Y_centroid', 'imageid', 'cell_id']
    marker_cols = list(marker_mapping.values())
    metadata_cols = ['Area', 'Classification']  # Keep some metadata

    # Check which columns actually exist
    available_coords = [col for col in coordinate_cols if col in df.columns]
    available_markers = [col for col in marker_cols if col in df.columns]
    available_metadata = [col for col in metadata_cols if col in df.columns]

    selected_cols = available_coords + available_markers + available_metadata
    df_clean = df[selected_cols].copy()

    # Step 4: Handle missing values
    print(f"📊 Data shape after cleaning: {df_clean.shape}")
    print(f"📊 Available markers: {available_markers}")

    # Check for missing values
    missing_counts = df_clean[available_markers].isnull().sum()
    if missing_counts.sum() > 0:
        print("⚠️  Missing values detected:")
        for marker, count in missing_counts.items():
            if count > 0:
                print(f"   {marker}: {count} missing values")

        # Fill missing values with 0 (assuming missing = no expression)
        df_clean[available_markers] = df_clean[available_markers].fillna(0)

    # Step 5: Create AnnData object
    # Expression data (markers only)
    X = df_clean[available_markers].values

    # Observation metadata (cell-level info)
    obs = df_clean[available_coords + available_metadata].copy()

    # Variable metadata (marker-level info)
    var = pd.DataFrame(index=available_markers)
    var['marker_type'] = 'protein'

    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs, var=var)

    # Step 6: Basic data validation
    print("\n🔍 Data validation:")
    print(f"   Total cells: {adata.n_obs:,}")
    print(f"   Total markers: {adata.n_vars}")
    print(f"   Images/regions: {adata.obs['imageid'].nunique()}")

    # Check coordinate ranges
    x_range = adata.obs['X_centroid'].max() - adata.obs['X_centroid'].min()
    y_range = adata.obs['Y_centroid'].max() - adata.obs['Y_centroid'].min()
    print(f"   Coordinate ranges: X={x_range:.0f}, Y={y_range:.0f}")

    # Check expression ranges
    expr_stats = pd.DataFrame(adata.X).describe()
    print(f"   Expression range: {expr_stats.loc['min'].min():.2f} - {expr_stats.loc['max'].max():.2f}")

    return adata


# =============================================================================
# 3. CREATE UPDATED PHENOTYPE WORKFLOW
# =============================================================================

def create_updated_phenotype_workflow():
    """
    Create phenotype_workflow.csv with your actual marker names
    """

    # Corrected phenotype workflow data
    phenotype_data = {
        'phenotype': [
            'M1_macrophages',
            'M2_macrophages',
            'Helper_T_cells',
            'Cytotoxic_T_cells',
            'Exhausted_Helper_T_cells',
            'Exhausted_Cytotoxic_T_cells',
            'Tumor_cells',
            'Tumor_proliferating_cells',
            'Antigen_presenting_cells_CD11c',
            'Antigen_presenting_cells_CD68_CD163_pos',
            'Antigen_presenting_cells_CD68_CD163_neg',
            'B_cells',
            'Vessels',
            'Tregs',
            'Other_immune_cells'
        ],
        'CD3': ['', '', 'pos', 'pos', 'pos', 'pos', '', '', '', '', '', '', '', 'pos', 'pos'],
        'CD4': ['', '', 'pos', '', 'pos', '', '', '', '', '', '', '', '', 'pos', ''],
        'CD8': ['', '', '', 'pos', '', 'pos', '', '', '', '', '', '', '', '', ''],
        'CD68': ['', '', '', '', '', '', '', '', '', 'pos', 'pos', '', '', '', ''],
        'CD163': ['pos', 'pos', '', '', '', '', '', '', '', 'pos', 'neg', '', '', '', ''],
        'PD1': ['neg', 'pos', '', '', 'pos', 'pos', '', '', '', '', '', '', '', '', ''],
        'PanCK': ['', '', '', '', '', '', '', 'pos', '', 'pos', 'pos', '', '', '', ''],
        'Ki67': ['', '', '', '', '', '', 'pos', 'pos', '', '', '', '', '', '', ''],
        'MHC-II': ['', '', '', '', '', '', '', '', 'pos', '', '', '', '', '', ''],
        'CD11c': ['', '', '', '', '', '', '', '', 'pos', '', '', '', '', '', ''],
        'CD20': ['', '', '', '', '', '', '', '', '', '', '', 'pos', '', '', ''],
        'CD31': ['', '', '', '', '', '', '', '', '', '', '', '', 'pos', '', ''],
        'FOXP3': ['', '', '', '', '', '', '', '', '', '', '', '', '', 'pos', '']
    }

    # Create DataFrame
    phenotype_df = pd.DataFrame(phenotype_data)

    # Add the first column as integers (required by scimap)
    phenotype_df.insert(0, 'id', range(1, len(phenotype_df) + 1))

    return phenotype_df


# =============================================================================
# 4. QUALITY CONTROL AND PREPROCESSING
# =============================================================================

def perform_quality_control(adata, slug="undefined"):
    """
    Perform comprehensive quality control checks
    """
    print("\n🔍 QUALITY CONTROL ANALYSIS")
    print("=" * 50)

    # Basic statistics
    print("📈 Basic Statistics:")
    print(f"Total cells: {adata.n_obs:,}")
    print(f"Total markers: {adata.n_vars}")
    print(f"Images/regions: {adata.obs['imageid'].nunique()}")

    # Cells per image/region
    cells_per_image = adata.obs['imageid'].value_counts()
    print(
        f"Cells per region - Mean: {cells_per_image.mean():.0f}, Range: {cells_per_image.min()}-{cells_per_image.max()}")

    # Expression statistics
    expr_df = pd.DataFrame(adata.X, columns=adata.var_names)
    print(f"\n📊 Expression Statistics:")
    print(expr_df.describe())

    # Create comprehensive QC plots
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))

    # 1. Overall expression distribution
    expr_flat = adata.X.flatten()
    axes[0, 0].hist(expr_flat, bins=50, alpha=0.7, edgecolor='black')
    axes[0, 0].set_title('Overall Expression Distribution')
    axes[0, 0].set_xlabel('Expression Level')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_yscale('log')

    # 2. Cells per region
    cells_per_image.plot(kind='bar', ax=axes[0, 1])
    axes[0, 1].set_title('Cells per Region')
    axes[0, 1].set_xlabel('Region ID')
    axes[0, 1].set_ylabel('Number of Cells')
    axes[0, 1].tick_params(axis='x', rotation=45)

    # 3. Marker expression boxplot
    expr_df.boxplot(ax=axes[0, 2])
    axes[0, 2].set_title('Marker Expression Distribution')
    axes[0, 2].tick_params(axis='x', rotation=45)
    axes[0, 2].set_yscale('log')

    # 4. Correlation heatmap
    corr_matrix = expr_df.corr()
    im1 = axes[1, 0].imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    axes[1, 0].set_title('Marker Correlation Matrix')
    axes[1, 0].set_xticks(range(len(adata.var_names)))
    axes[1, 0].set_yticks(range(len(adata.var_names)))
    axes[1, 0].set_xticklabels(adata.var_names, rotation=45)
    axes[1, 0].set_yticklabels(adata.var_names)
    plt.colorbar(im1, ax=axes[1, 0])

    # 5. Spatial distribution
    scatter = axes[1, 1].scatter(adata.obs['X_centroid'], adata.obs['Y_centroid'],
                                 s=0.5, alpha=0.6, c='blue')
    axes[1, 1].set_title('Spatial Distribution of Cells')
    axes[1, 1].set_xlabel('X Coordinate')
    axes[1, 1].set_ylabel('Y Coordinate')
    axes[1, 1].set_aspect('equal')

    # 6. Cell area distribution (if available)
    if 'Area' in adata.obs.columns:
        adata.obs['Area'].hist(bins=50, ax=axes[1, 2])
        axes[1, 2].set_title('Cell Area Distribution')
        axes[1, 2].set_xlabel('Area (pixels²)')
        axes[1, 2].set_ylabel('Frequency')
    else:
        axes[1, 2].text(0.5, 0.5, 'Area data\nnot available',
                        ha='center', va='center', transform=axes[1, 2].transAxes)
        axes[1, 2].set_title('Cell Area Distribution')

    # 7-9. Individual marker distributions for key markers
    key_markers = ['CD3', 'CD68', 'PanCK']
    for i, marker in enumerate(key_markers):
        if marker in adata.var_names:
            marker_idx = adata.var_names.get_loc(marker)
            axes[2, i].hist(adata.X[:, marker_idx], bins=50, alpha=0.7)
            axes[2, i].set_title(f'{marker} Expression')
            axes[2, i].set_xlabel('Expression Level')
            axes[2, i].set_ylabel('Frequency')
        else:
            axes[2, i].text(0.5, 0.5, f'{marker}\nnot available',
                            ha='center', va='center', transform=axes[2, i].transAxes)
            axes[2, i].set_title(f'{marker} Expression')

    plt.tight_layout()
    plt.savefig(f'normalized_qc/quality_control_plots_{slug}.png', dpi=300, bbox_inches='tight')
    #plt.show()

    # Summary statistics
    print("\n📋 QC Summary:")
    zero_expr_cells = (adata.X.sum(axis=1) == 0).sum()
    print(f"   Cells with zero expression: {zero_expr_cells} ({zero_expr_cells / adata.n_obs * 100:.1f}%)")

    high_expr_cells = (adata.X.sum(axis=1) > np.percentile(adata.X.sum(axis=1), 95)).sum()
    print(
        f"   High expression cells (>95th percentile): {high_expr_cells} ({high_expr_cells / adata.n_obs * 100:.1f}%)")

    return adata


# =============================================================================
# 5. GATE IDENTIFICATION AND DATA RESCALING
# =============================================================================

def perform_gating_and_rescaling(adata, method='auto'):
    """
    Perform gating and data rescaling

    Parameters:
    -----------
    method : str
        'auto' for automated gating, 'manual' for manual gates
    """
    print(f"\n⚖️  GATING AND RESCALING ({method.upper()})")
    print("=" * 50)

    if method == 'manual':
        print("🎯 For manual gating, you would typically:")
        print("1. Use sm.pl.gate_finder() for each marker")
        print("2. Visually determine optimal gate values")
        print("3. Create manual_gates.csv file")
        print("4. For now, using automated gating...")

    try:
        # Automated gating using Gaussian Mixture Models
        print("🤖 Performing automated gating using GMM...")
        adata = sm.pp.rescale(adata)

        print("✅ Data rescaling completed")

        # Visualize rescaling results
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))

        # Plot before/after distributions for key markers
        key_markers = ['CD3', 'CD68', 'PanCK', 'CD4', 'CD8', 'Ki67']
        available_key_markers = [m for m in key_markers if m in adata.var_names][:6]

        for i, marker in enumerate(available_key_markers):
            row = i // 3
            col = i % 3

            marker_idx = adata.var_names.get_loc(marker)
            axes[row, col].hist(adata.X[:, marker_idx], bins=50, alpha=0.7,
                                color='blue', label='After Rescaling')
            axes[row, col].set_title(f'{marker} - After Rescaling')
            axes[row, col].set_xlabel('Expression Level (Rescaled)')
            axes[row, col].set_ylabel('Frequency')
            axes[row, col].axvline(x=0.5, color='red', linestyle='--',
                                   label='Default Gate (0.5)')
            axes[row, col].legend()

        # Hide empty subplots
        for i in range(len(available_key_markers), 6):
            row = i // 3
            col = i % 3
            axes[row, col].set_visible(False)

        # extract basename from data_path
        import os
        base_name = os.path.basename(args.data_path).replace('.csv', '')
        slug = base_name.replace('cell_data_wih_regions_', '')

        plt.tight_layout()
        plt.savefig(f'normalized_qc/rescaling_results_{slug}.png', dpi=300, bbox_inches='tight')
        # plt.show()

        # Print rescaling summary
        print(f"\n📊 Rescaling Summary:")
        rescaled_stats = pd.DataFrame(adata.X, columns=adata.var_names).describe()
        print(
            f"   Expression range after rescaling: {rescaled_stats.loc['min'].min():.3f} - {rescaled_stats.loc['max'].max():.3f}")
        print(f"   Mean expression: {rescaled_stats.loc['mean'].mean():.3f}")

    except Exception as e:
        print(f"❌ Error during rescaling: {e}")
        print("💡 This might be due to data format issues. Check your marker columns.")
        raise e

    return adata


# =============================================================================
# 6. CELL PHENOTYPING
# =============================================================================

def perform_phenotyping(adata):
    """
    Perform cell phenotyping using the workflow CSV
    """
    print("\n🏷️  CELL PHENOTYPING")
    print("=" * 50)

    try:
        # Use the corrected phenotype workflow
        phenotype_path = 'data/phenotype_workflow_corrected.csv'

        # Create corrected phenotype workflow if it doesn't exist
        #if not os.path.exists(phenotype_path):
        #    print("Creating corrected phenotype workflow...")
        #    phenotype = create_updated_phenotype_workflow()
        #    phenotype.to_csv(phenotype_path, index=False)
        #else:
        phenotype = pd.read_csv(phenotype_path)

        #print(f"📋 Loaded phenotype workflow with {len(phenotype)} phenotypes")

        # Ensure proper data types
        #phenotype = phenotype.fillna('')  # Replace NaN with empty strings

        # Make sure first column is numeric
        #if 'unnamed' in phenotype.columns:
        #    phenotype = phenotype.drop('unnamed', axis=1)

        # Ensure we have a proper numeric first column
        #if not phenotype.columns[0] == 'id':
        #    phenotype.insert(0, 'id', range(1, len(phenotype) + 1))

        # Convert all marker columns to strings
        marker_columns = phenotype.columns[2:]  # Skip 'id' and 'phenotype'
        for col in marker_columns:
            phenotype[col] = phenotype[col].astype(str).replace('nan', '')

        # Display workflow for verification
        print("\n📋 Phenotype workflow preview:")
        print(phenotype.head())

        # Perform phenotyping
        adata = sm.tl.phenotype_cells(adata,
                                      phenotype=phenotype,
                                      label="phenotype")
        # Analyze phenotyping results
        phenotype_counts = adata.obs['phenotype'].value_counts()
        print("\n📊 Phenotyping Results:")
        print(phenotype_counts)

        # Calculate success metrics
        total_cells = len(adata.obs)
        unknown_cells = phenotype_counts.get('Unknown', 0)
        success_rate = (total_cells - unknown_cells) / total_cells * 100

        print(f"\n🎯 Phenotyping Success Rate: {success_rate:.1f}%")
        print(f"   Successfully classified: {total_cells - unknown_cells:,} cells")
        print(f"   Unknown cells: {unknown_cells:,} cells ({unknown_cells / total_cells * 100:.1f}%)")

        # Create comprehensive visualization
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))

        # 1. Pie chart of phenotypes
        # Only show top phenotypes to avoid clutter
        top_phenotypes = phenotype_counts.head(10)
        if len(phenotype_counts) > 10:
            other_count = phenotype_counts.iloc[10:].sum()
            top_phenotypes['Others'] = other_count

        wedges, texts, autotexts = axes[0, 0].pie(top_phenotypes.values,
                                                  labels=top_phenotypes.index,
                                                  autopct='%1.1f%%',
                                                  startangle=90)
        axes[0, 0].set_title('Cell Phenotype Distribution (Top 10)')

        # 2. Bar chart
        phenotype_counts.head(15).plot(kind='bar', ax=axes[0, 1])
        axes[0, 1].set_title('Cell Counts by Phenotype (Top 15)')
        axes[0, 1].tick_params(axis='x', rotation=45)
        axes[0, 1].set_ylabel('Number of Cells')

        # 3. Phenotype by region
        if 'imageid' in adata.obs.columns and adata.obs['imageid'].nunique() > 1:
            phenotype_by_region = pd.crosstab(adata.obs['imageid'],
                                              adata.obs['phenotype'])
            # Show only top phenotypes
            top_pheno_names = phenotype_counts.head(8).index
            phenotype_by_region[top_pheno_names].plot(kind='bar', stacked=True, ax=axes[0, 2])
            axes[0, 2].set_title('Phenotype Distribution by Region')
            axes[0, 2].tick_params(axis='x', rotation=45)
            axes[0, 2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            axes[0, 2].text(0.5, 0.5, 'Single region\nor no region data',
                            ha='center', va='center', transform=axes[0, 2].transAxes)
            axes[0, 2].set_title('Phenotype Distribution by Region')

        # 4. Classification success
        success_data = {'Successfully\nClassified': success_rate, 'Unknown': 100 - success_rate}
        bars = axes[1, 0].bar(success_data.keys(), success_data.values(),
                              color=['lightgreen', 'orange'])
        axes[1, 0].set_title(f'Classification Success\n({success_rate:.1f}% Success Rate)')
        axes[1, 0].set_ylabel('Percentage of Cells')
        axes[1, 0].set_ylim(0, 100)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            axes[1, 0].annotate(f'{height:.1f}%',
                                xy=(bar.get_x() + bar.get_width() / 2, height),
                                xytext=(0, 3),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='bottom')

        # 5. Immune vs non-immune breakdown
        immune_phenotypes = ['Helper_T_cells', 'Cytotoxic_T_cells', 'Exhausted_Helper_T_cells',
                             'Exhausted_Cytotoxic_T_cells', 'M1_macrophages', 'M2_macrophages',
                             'B_cells', 'Tregs', 'Other_immune_cells', 'Antigen_presenting_cells_CD11c',
                             'Antigen_presenting_cells_CD68_CD163_pos', 'Antigen_presenting_cells_CD68_CD163_neg']

        immune_count = adata.obs[adata.obs['phenotype'].isin(immune_phenotypes)].shape[0]
        tumor_count = adata.obs[adata.obs['phenotype'].str.contains('Tumor', na=False)].shape[0]
        vessel_count = adata.obs[adata.obs['phenotype'] == 'Vessels'].shape[0]
        unknown_count = adata.obs[adata.obs['phenotype'] == 'Unknown'].shape[0]
        other_count = total_cells - immune_count - tumor_count - vessel_count - unknown_count

        categories = ['Immune', 'Tumor', 'Vessels', 'Unknown', 'Other']
        counts = [immune_count, tumor_count, vessel_count, unknown_count, other_count]

        axes[1, 1].pie(counts, labels=categories, autopct='%1.1f%%', startangle=90)
        axes[1, 1].set_title('Major Cell Category Distribution')

        # 6. Expression heatmap for top phenotypes
        if len(phenotype_counts) > 1:
            top_phenotypes_for_heatmap = phenotype_counts.head(8).index
            phenotype_expr = []
            phenotype_names = []

            for pheno in top_phenotypes_for_heatmap:
                if pheno != 'Unknown':  # Skip unknown for cleaner visualization
                    cells_of_type = adata.obs['phenotype'] == pheno
                    if cells_of_type.sum() > 0:
                        mean_expr = adata.X[cells_of_type].mean(axis=0)
                        phenotype_expr.append(mean_expr)
                        phenotype_names.append(pheno)

            if phenotype_expr:
                expr_matrix = np.array(phenotype_expr)
                im = axes[1, 2].imshow(expr_matrix, aspect='auto', cmap='viridis')
                axes[1, 2].set_title('Mean Expression by Phenotype')
                axes[1, 2].set_yticks(range(len(phenotype_names)))
                axes[1, 2].set_yticklabels(phenotype_names)
                axes[1, 2].set_xticks(range(len(adata.var_names)))
                axes[1, 2].set_xticklabels(adata.var_names, rotation=45)
                plt.colorbar(im, ax=axes[1, 2])

        plt.tight_layout()
        plt.savefig('normalized_qc/phenotyping_results.png', dpi=300, bbox_inches='tight')
        #plt.show()

        # Save detailed results
        phenotype_summary = pd.DataFrame({
            'Phenotype': phenotype_counts.index,
            'Count': phenotype_counts.values,
            'Percentage': (phenotype_counts.values / total_cells * 100).round(2)
        })
        phenotype_summary.to_csv('phenotype_summary.csv', index=False)
        print(f"💾 Detailed results saved to phenotype_summary.csv")

    except Exception as e:
        print(f"❌ Error during phenotyping: {e}")
        print("💡 Check your phenotype_workflow.csv file and marker names")
        raise e

    return adata


# =============================================================================
# 7. SPATIAL ANALYSIS
# =============================================================================

def perform_spatial_analysis(adata):
    """
    Comprehensive spatial analysis adapted for your data
    """
    print("\n🗺️  SPATIAL ANALYSIS")
    print("=" * 50)

    # Check coordinate scale and provide guidance
    x_range = adata.obs['X_centroid'].max() - adata.obs['X_centroid'].min()
    y_range = adata.obs['Y_centroid'].max() - adata.obs['Y_centroid'].min()
    avg_range = (x_range + y_range) / 2

    # Suggest appropriate radius based on coordinate scale
    suggested_radius = max(50, avg_range / 100)  # Adaptive radius
    suggested_k = min(15, max(5, int(adata.n_obs / 1000)))  # Adaptive k

    print(f"📏 Coordinate analysis:")
    print(f"   Image dimensions: {x_range:.0f} x {y_range:.0f}")
    print(f"   Suggested radius for spatial analysis: {suggested_radius:.0f}")
    print(f"   Suggested k-neighbors: {suggested_k}")

    try:
        # 7.1 Spatial distribution visualization
        print("\n📍 Creating spatial distribution plots...")

        # Create spatial plots
        sm.pl.spatial_scatterPlot(adata,
                                  colorBy=['phenotype'],
                                  s=2,  # Smaller points for dense data
                                  figsize=(12, 10),
                                  saveDir='./plots/',
                                  fileName='spatial_phenotype_map.pdf')

        # Spatial plots for key markers
        key_markers = ['CD3', 'CD68', 'PanCK', 'Ki67', 'CD4', 'CD8']
        available_markers = [m for m in key_markers if m in adata.var_names]

        if available_markers:
            sm.pl.spatial_scatterPlot(adata,
                                      colorBy=available_markers[:4],  # Limit to 4 for clarity
                                      s=2,
                                      figsize=(15, 12),
                                      saveDir='./plots/',
                                      fileName='spatial_markers_map.pdf')

        # 7.2 Neighborhood analysis
        print("\n🏘️  Performing neighborhood analysis...")

        # Calculate spatial neighbors
        adata = sm.tl.spatial_knn(adata,
                                  x_coordinate='X_centroid',
                                  y_coordinate='Y_centroid',
                                  k=suggested_k,
                                  method='ball_tree')

        # Spatial lag analysis
        adata = sm.tl.spatial_lda(adata,
                                  phenotype='phenotype',
                                  k=suggested_k)

        # 7.3 Spatial interaction analysis
        print("\n🤝 Analyzing spatial interactions...")

        adata = sm.tl.spatial_interaction(adata,
                                          phenotype='phenotype',
                                          method='radius',
                                          radius=suggested_radius)

        # 7.4 Proximity analysis for key cell types
        print("\n📏 Analyzing cell-cell proximity...")

        # Focus on clinically relevant interactions
        target_phenotypes = ['Tumor_cells', 'Helper_T_cells', 'Cytotoxic_T_cells', 'M1_macrophages']
        available_phenotypes = [p for p in target_phenotypes
                                if p in adata.obs['phenotype'].unique()]

        for phenotype in available_phenotypes:
            if adata.obs[adata.obs['phenotype'] == phenotype].shape[0] > 10:  # Only if enough cells
                adata = sm.tl.spatial_distance(adata,
                                               x_coordinate='X_centroid',
                                               y_coordinate='Y_centroid',
                                               phenotype='phenotype',
                                               subset=phenotype,
                                               label=f'distance_to_{phenotype}')

        # 7.5 Spatial statistics visualization
        print("\n📊 Creating spatial analysis visualizations...")

        # Spatial interaction heatmap
        if 'spatial_interaction' in adata.uns:
            interaction_matrix = adata.uns['spatial_interaction']

            plt.figure(figsize=(12, 10))
            sns.heatmap(interaction_matrix, annot=True, cmap='coolwarm', center=0, fmt='.2f')
            plt.title('Spatial Interaction Analysis\n(Positive = Attraction, Negative = Avoidance)')
            plt.xlabel('Cell Type')
            plt.ylabel('Cell Type')
            plt.tight_layout()
            plt.savefig('normalized_qc/spatial_interactions.png', dpi=300, bbox_inches='tight')
            #plt.show()

            # Save interaction results
            interaction_matrix.to_csv('spatial_interactions.csv')

        print("✅ Spatial analysis completed successfully!")

    except Exception as e:
        print(f"⚠️  Spatial analysis error: {e}")
        print(f"💡 Try adjusting radius (current: {suggested_radius}) or k (current: {suggested_k}) parameters")
        print("💡 Or check if you have sufficient cells for spatial analysis")

    return adata

# =============================================================================
# 8. STATISTICAL ANALYSIS
# =============================================================================

def perform_statistical_analysis(adata):
    """
    Statistical analysis and hypothesis testing
    """
    print("\n📊 STATISTICAL ANALYSIS")
    print("=" * 50)

    try:
        # 8.1 Phenotype abundance analysis
        print("🧮 Analyzing phenotype abundances...")

        if adata.obs['imageid'].nunique() > 1:
            # Calculate phenotype proportions per region
            phenotype_props = adata.obs.groupby(['imageid', 'phenotype']).size().unstack(fill_value=0)

            # Convert to proportions
            phenotype_props_pct = phenotype_props.div(phenotype_props.sum(axis=1), axis=0)

            print("📈 Phenotype proportion statistics:")
            print(phenotype_props_pct.describe())

            # Save results
            phenotype_props.to_csv('phenotype_counts_by_region.csv')
            phenotype_props_pct.to_csv('phenotype_proportions_by_region.csv')

            # Statistical comparisons (if you have treatment groups)
            print("\n🔬 Statistical comparisons between regions:")
            for phenotype in phenotype_props_pct.columns:
                if phenotype != 'Unknown':
                    values = phenotype_props_pct[phenotype].dropna()
                    if len(values) > 1:
                        mean_prop = values.mean()
                        std_prop = values.std()
                        print(f"   {phenotype}: {mean_prop:.3f} ± {std_prop:.3f}")

            # Visualization
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))

            # Heatmap of proportions
            sns.heatmap(phenotype_props_pct.T, annot=True, fmt='.3f', ax=axes[0, 0], cmap='viridis')
            axes[0, 0].set_title('Phenotype Proportions by Region')
            axes[0, 0].set_xlabel('Region')
            axes[0, 0].set_ylabel('Phenotype')

            # Box plots of key phenotypes
            key_phenotypes = ['Helper_T_cells', 'Cytotoxic_T_cells', 'M1_macrophages', 'Tumor_cells']
            available_key = [p for p in key_phenotypes if p in phenotype_props_pct.columns]

            if available_key:
                phenotype_props_pct[available_key].boxplot(ax=axes[0, 1])
                axes[0, 1].set_title('Distribution of Key Phenotype Proportions')
                axes[0, 1].tick_params(axis='x', rotation=45)
                axes[0, 1].set_ylabel('Proportion')

            # Regional variation
            regional_cv = phenotype_props_pct.std() / phenotype_props_pct.mean()
            regional_cv.dropna().sort_values(ascending=False).head(10).plot(kind='bar', ax=axes[1, 0])
            axes[1, 0].set_title('Phenotype Variability Across Regions\n(Coefficient of Variation)')
            axes[1, 0].set_ylabel('CV (std/mean)')
            axes[1, 0].tick_params(axis='x', rotation=45)

        else:
            print("📊 Single region analysis:")
            phenotype_counts = adata.obs['phenotype'].value_counts()
            phenotype_props = phenotype_counts / phenotype_counts.sum()

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Single region visualization
            phenotype_counts.head(10).plot(kind='bar', ax=axes[0])
            axes[0].set_title('Phenotype Counts')
            axes[0].tick_params(axis='x', rotation=45)

        # 8.2 Spatial enrichment analysis
        print("\n🗺️  Spatial enrichment analysis...")

        if 'spatial_lda' in adata.obsm.keys():
            # Spatial enrichment scoring
            adata = sm.tl.spatial_pscore(adata,
                                         phenotype='phenotype',
                                         method='radius',
                                         radius=max(50, (x_range + y_range) / 200))

            if 'spatial_pscore' in adata.uns:
                spatial_results = adata.uns['spatial_pscore']
                print("📍 Spatial enrichment results (top 10):")
                print(spatial_results.head(10))

                # Save spatial results
                spatial_results.to_csv('spatial_enrichment_results.csv')

        # 8.3 Distance analysis
        print("\n📏 Distance-based analysis...")

        distance_cols = [col for col in adata.obs.columns if col.startswith('distance_to_')]

        if distance_cols:
            distance_df = adata.obs[distance_cols + ['phenotype']].copy()

            # Summary statistics
            print("📏 Distance summary statistics:")
            print(distance_df[distance_cols].describe())

            # Distance distributions by phenotype
            if len(distance_cols) > 0:
                fig, axes = plt.subplots(1, len(distance_cols), figsize=(6 * len(distance_cols), 5))
                if len(distance_cols) == 1:
                    axes = [axes]

                for i, dist_col in enumerate(distance_cols):
                    target_cell_type = dist_col.replace('distance_to_', '')

                    # Plot distance distributions for different phenotypes
                    phenotypes_to_plot = adata.obs['phenotype'].value_counts().head(5).index

                    for phenotype in phenotypes_to_plot:
                        if phenotype != target_cell_type:  # Don't plot distance to self
                            subset_distances = distance_df[distance_df['phenotype'] == phenotype][dist_col]
                            if len(subset_distances) > 10:
                                axes[i].hist(subset_distances, alpha=0.6, label=phenotype, bins=30)

                    axes[i].set_title(f'Distance to {target_cell_type}')
                    axes[i].set_xlabel('Distance (pixels)')
                    axes[i].set_ylabel('Frequency')
                    axes[i].legend()

                plt.tight_layout()
                plt.savefig('normalized_qc/distance_analysis.png', dpi=300, bbox_inches='tight')
                plt.show()

            # Save distance analysis
            distance_df.to_csv('distance_analysis_results.csv', index=False)

        # 8.4 Marker correlation analysis
        print("\n🔗 Marker correlation analysis...")

        marker_corr = np.corrcoef(adata.X.T)
        marker_corr_df = pd.DataFrame(marker_corr,
                                      index=adata.var_names,
                                      columns=adata.var_names)

        # Find strongest correlations
        corr_pairs = []
        for i in range(len(adata.var_names)):
            for j in range(i + 1, len(adata.var_names)):
                corr_pairs.append({
                    'Marker1': adata.var_names[i],
                    'Marker2': adata.var_names[j],
                    'Correlation': marker_corr[i, j]
                })

        corr_pairs_df = pd.DataFrame(corr_pairs)
        corr_pairs_df = corr_pairs_df.sort_values('Correlation', key=abs, ascending=False)

        print("🔗 Strongest marker correlations:")
        print(corr_pairs_df.head(10))

        # Visualize correlations
        plt.figure(figsize=(10, 8))
        mask = np.triu(np.ones_like(marker_corr, dtype=bool))
        sns.heatmap(marker_corr_df, mask=mask, annot=True, cmap='coolwarm', center=0,
                    square=True, fmt='.2f')
        plt.title('Marker Expression Correlations')
        plt.tight_layout()
        plt.savefig('normalized_qc/marker_correlations.png', dpi=300, bbox_inches='tight')
        plt.show()

        # Save correlation results
        marker_corr_df.to_csv('marker_correlations.csv')
        corr_pairs_df.to_csv('marker_correlation_pairs.csv', index=False)

        print("✅ Statistical analysis completed!")

    except Exception as e:
        print(f"⚠️  Statistical analysis error: {e}")
        print("💡 Check data format and ensure sufficient sample size")

    return adata


# =============================================================================
# 9. MAIN WORKFLOW FUNCTION
# =============================================================================

def run_complete_analysis(data_path):
    """
    Run the complete SCIMAP analysis workflow adapted for your data format

    Parameters:
    -----------
    data_path : str
        Path to your CSV file with QuPath data
    """

    print("🚀 STARTING COMPLETE SCIMAP ANALYSIS")
    print("🔬 Adapted for QuPath-exported IHC data")
    print("=" * 60)

    try:
        # Step 1: Load and preprocess data
        adata = load_and_preprocess_data(data_path)

        # Step 2: Create/update phenotype workflow
        #create_updated_phenotype_workflow()

        # Step 3: Quality control
        slug = data_path.split('/')[-1].replace('.csv', '')
        adata = perform_quality_control(adata, slug)

        return adata

    except Exception as e:
        print(f"\n❌ ANALYSIS FAILED: {e}")
        print("💡 Check the troubleshooting section in the script")
        import traceback
        traceback.print_exc()
        return None


import argparse
# Load sample name from command
parser = argparse.ArgumentParser(description="Run complete phenotyping workflow")
parser.add_argument('--data-path', type=str, default=None,
                  help='Path to the cell data CSV file')
args = parser.parse_args()

if __name__ == '__main__':
    if args.data_path is None:
        data_path_list = glob("cell_data/normalized/*.csv")
        for data_path in data_path_list:
            adata = run_complete_analysis(data_path)
    else:
        adata = run_complete_analysis(args.data_path)
    print()