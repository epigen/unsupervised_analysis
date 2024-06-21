#### generate diagnostic plots for Uniform Manifold Approximation and Projection (UMAP) ####

#### libraries
# general
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
# dimensionality reduction
import umap
import umap.plot

#### configurations

# ipnuts
object_path = snakemake.input["umap_object"]

# outputs
plot_connectivity_path = snakemake.output["connectivity_plot"]

### load data
with open(object_path, 'rb') as f:
    umap_obj = pickle.load(f)
    
### generate & save UMAP connectivity plot

# umap.plot.connectivity(umap_obj, show_points=True)
umap.plot.connectivity(umap_obj, edge_bundling='hammer').figure.savefig(plot_connectivity_path) 
