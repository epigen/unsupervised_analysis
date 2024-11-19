#### determine distance matrix using scipy ####

#### libraries
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import math
# from fastdist import fastdist

#### configurations

# ipnuts
data_path = snakemake.input["data"]

# outputs
distance_matrix_path = snakemake.output["distance_matrix"]

# parameters
samples_by_features = int(snakemake.params['samples_by_features'])
data_or_feature = snakemake.wildcards["type"] # "observations" or "features"
metric = snakemake.wildcards["metric"]
n_observations = snakemake.config["heatmap"]["n_observations"]
n_features = snakemake.config["heatmap"]["n_features"]

### load data

# check data orientation to fit: samples/observations x features
data = pd.read_csv(data_path, index_col=0)
if (samples_by_features == 0 and data_or_feature == "observations") or (samples_by_features == 1 and data_or_feature == "features"):
    data = data.T


# retain only highly variable features
if data_or_feature == "features":
    variances = data.var(axis=1)
    if isinstance(n_features, float) or n_features==1:
        n_features = int(math.floor(n_features * data.shape[0]))
    top_features = variances.nlargest(n_features).index
    data = data.loc[top_features,:]

# downsample observations
if data_or_feature == "observations":
    if isinstance(n_observations, float) or n_observations==1:
        n_observations = int(math.floor(n_observations * data.shape[0]))
    if n_observations < data.shape[0]:
        data = data.sample(n=n_observations, random_state=42)
    
# Convert DataFrame to NumPy array
# data_np = data.to_numpy()

### Distance matrix calculation
# Pairwise distances between observations in n-dimensional space.

# scipy
dist_mtx = pdist(data, metric=metric)

# fastdist
# metric_function = getattr(fastdist, metric)
# dist_mtx = fastdist.matrix_pairwise_distance(data_np, metric_function, metric, return_matrix=True) 

# convert to squareform dataframe
dist_mtx_df = pd.DataFrame(squareform(dist_mtx), index=data.index, columns=data.index)
# dist_mtx_df = pd.DataFrame(dist_mtx, index=data.index, columns=data.index)

### save distance matrix

dist_mtx_df.to_csv(distance_matrix_path)
