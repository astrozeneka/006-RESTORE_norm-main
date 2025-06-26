

import os

import scimap as sm
import anndata as ad
from matplotlib import pyplot as plt
import scimap as sm
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path


if __name__ == '__main__':
    os.mkdir("results-phenotypes/distance_to_pheno") if not os.path.exists("results-phenotypes/distance_to_pheno") else None
    data_path = "results-phenotypes/phenotypes/S19_12126B1.h5ad"

    # Read the data
    adata = ad.read_h5ad(data_path)
    adata.obs['imageid'] = 'image_1'  # hard-code the image-id

    # Compute the spatial interaction  using spatial_interaction
    si_adata = sm.tl.spatial_interaction(
        adata,
        x_coordinate='X_centroid',
        y_coordinate='Y_centroid',
        phenotype='phenotype',
        method='radius',
        radius=250, # need to adapt
        permutation=100, # Need to adapt
        pval_method='zscore' # default is histocat
    )

    # Basic interaction heatmap (doesn't work)
    #plt.figure(figsize=(12, 10))
    # Neighborhood interaction heatmap
    sm.pl.spatial_interaction(
        si_adata,
        spatial_interaction='spatial_interaction',
        summarize_plot=True, # only "True" for single image
        p_val=0.05,
        row_cluster=False,
        col_cluster=False,
        cmap='viridis',
        nonsig_color='white',
        fileName="spatial_interaction.png",
        saveDir="results-phenotypes/distance_to_pheno"
    )
    #plt.tight_layout()

    print()

    print()