

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
    data_path = "results-phenotypes/phenotypes/S19_12126B1.h5ad"

    # Read the data
    adata = ad.read_h5ad(data_path)
    adata.obs['imageid'] = 'image_1'  # hard-code the image-id

    # Compute the distance to phenotype using spatial_distance
    distance_adata = sm.tl.spatial_distance(
        adata,
        x_coordinate='X_centroid',
        y_coordinate='Y_centroid',
        phenotype='phenotype'    )

    print()

    # summary heatmap
    sm.pl.spatial_distance(
        adata,
        heatmap_cmap='viridis'
    )

    sm.pl.spatial_distance(adata, method='distribution', distance_from='Exhausted_Helper_T_cells', imageid='ImageId',
                           x_axis="distance", y_axis="imageid", plot_type="kde")

    print()