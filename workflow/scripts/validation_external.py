#### load all clustering results and metadata and determine external cluster indices ####

#### libraries
import os
import numpy as np
import pandas as pd
from sklearn import metrics

#### configurations

# inputs
metadata_path = os.path.join(snakemake.input["metadata"])
clusterings_path = os.path.join(snakemake.input["clusterings"])

# outputs
result_paths = snakemake.output

# load the clustering results and categorical metadata
metadata = pd.read_csv(metadata_path, index_col=0)
metadata = metadata.fillna('')

clustering_results = pd.read_csv(clusterings_path, index_col=0)

indices = [s.split('external_index_')[1].split('.csv')[0] for s in result_paths]

# identify categorical metadata
meta_cat = []
for variable in metadata.columns:
    unique_vals = list(metadata[variable].unique())
    
    # check if integer AND less than 25 unique values -> categorical metadata
    if all([isinstance(i, (int, np.int64)) for i in unique_vals]) and len(unique_vals)<25:
        meta_cat.append(variable)
        continue
    
    if all([isinstance(i, (str, bool, np.bool_)) for i in unique_vals]):
        meta_cat.append(variable)
        
    else:
        print("variable {} not categorical".format(variable))

        # subset for categorical data
categorical_metadata = metadata.loc[:,meta_cat]

# Ensure that the clustering results and categorical metadata have the same indices
# fix metadata indices if they do not agree with data as they come from outside the workflow (e.g., R)
if not(set(clustering_results.index) == set(categorical_metadata.index)):
    if categorical_metadata.index.inferred_type=='string':
        categorical_metadata.index = [idx.replace('-','.') for idx in categorical_metadata.index]

# Reorder the rows of the clustering results DataFrame
clustering_results = clustering_results.reindex(categorical_metadata.index)

idx_dfs = {}

# Create a DataFrame for each index
for idx in indices:
    idx_dfs[idx] = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)
    
# For each clustering result
for clustering in clustering_results.columns:
    # For each categorical metadata
    for metadata in categorical_metadata.columns:
        # Calculate & store the scores
        idx_dfs["AMI"].loc[clustering, metadata] = metrics.adjusted_mutual_info_score(categorical_metadata[metadata], clustering_results[clustering])
        idx_dfs["ARI"].loc[clustering, metadata] = metrics.adjusted_rand_score(categorical_metadata[metadata], clustering_results[clustering])
        idx_dfs["FMI"].loc[clustering, metadata] = metrics.fowlkes_mallows_score(categorical_metadata[metadata], clustering_results[clustering])
        
        homogeneity, completeness, v_measure = metrics.homogeneity_completeness_v_measure(categorical_metadata[metadata], clustering_results[clustering])
        idx_dfs["Homogeneity"].loc[clustering, metadata] = homogeneity
        idx_dfs["Completeness"].loc[clustering, metadata] = completeness
        idx_dfs["V"].loc[clustering, metadata] = v_measure

        
# # Save the DataFrames as CSV files
for i, idx in enumerate(indices):
    idx_dfs[idx].to_csv(result_paths[i])
