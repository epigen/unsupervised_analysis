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
object_path = snakemake.input["umap_object"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP/UMAP_correlation_100_0.1_2_object.pickle"

# outputs
plot_connectivity_path = snakemake.output["connectivity_plot"] # "/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP/plots/UMAP_correlation_100_0.1_2_connectivity.png"

result_dir = os.path.dirname(plot_connectivity_path)

# make directory if not existing
if not os.path.exists(result_dir):
    os.makedirs(result_dir, exist_ok=True)

### load data
with open(object_path, 'rb') as f:
    umap_obj = pickle.load(f)
    
### generate & save UMAP connectivity plot

# umap.plot.connectivity(umap_obj, show_points=True)
umap.plot.connectivity(umap_obj, edge_bundling='hammer').figure.savefig(plot_connectivity_path) 
