#### perform knn graph generation from Uniform Manifold Approximation and Projection (UMAP) ####

#### libraries
# general
import os
import pickle
import pandas as pd
# dimensionality reduction
import umap
from umap.umap_ import nearest_neighbors

#### configurations

# ipnuts
data_path = snakemake.input[0] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/KOcall_NonTargeting/counts/CORRECTED_RNA.csv"
# outputs
result_object_path = snakemake.output["result_object"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP_correlation_100_object.pickle"

result_dir = os.path.dirname(result_object_path)

# parameters
samples_by_features = snakemake.params['samples_by_features'] #0
metric = snakemake.params['metric'] # "correlation"
n_neighbors = int(snakemake.params['n_neighbors']) #100

# make directory if not existing
if not os.path.exists(result_dir):
    os.makedirs(result_dir, exist_ok=True)

### load data

# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T
    
### get knn-graph
knn = nearest_neighbors(data,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        metric_kwds=None,
                        angular=False,
                        random_state=42,
                        low_memory=True, 
                        use_pynndescent=True, 
                        n_jobs=-1, 
                        verbose=False
                       )

## save knn graph object
with open(result_object_path, 'wb') as f:
    pickle.dump(knn, f, pickle.HIGHEST_PROTOCOL)