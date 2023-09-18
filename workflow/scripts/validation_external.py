#### load all clustering results and metadata and determine external cluster indices ####

#### libraries
import os
import numpy as np
import pandas as pd
from sklearn import metrics

#### configurations

# inputs
metadata_path = os.path.join(snakemake.input["metadata"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/data/digits_labels.csv"
clusterings_path = os.path.join(snakemake.input["clusterings"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/metadata_clusterings.csv"

# outputs
# ami_path = os.path.join(snakemake.output["ami"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_AMI.csv"
# ari_path = os.path.join(snakemake.output["ari"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_ARI.csv"
# fmi_path = os.path.join(snakemake.output["fmi"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_FMI.csv"
# homogeneity_path = os.path.join(snakemake.output["homogeneity"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_Homogeneity.csv"
# completeness_path = os.path.join(snakemake.output["completeness"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_Completeness.csv"
# v_path = os.path.join(snakemake.output["v"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/cluster_validation/external_index_V.csv"

result_paths = snakemake.output #  [".test/results/unsupervised_analysis/digits/cluster_validation/external_index_AMI.csv", ".test/results/unsupervised_analysis/digits/cluster_validation/external_index_ARI.csv", ".test/results/unsupervised_analysis/digits/cluster_validation/external_index_FMI.csv", ".test/results/unsupervised_analysis/digits/cluster_validation/external_index_Homogeneity.csv", ".test/results/unsupervised_analysis/digits/cluster_validation/external_index_Completeness.csv", ".test/results/unsupervised_analysis/digits/cluster_validation/external_index_V.csv"]


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

# # Create a DataFrame for each index
# ami_scores = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)
# ari_scores = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)
# fmi_scores = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)
# homogeneity_scores = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)
# completeness_scores = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)
# v_measure_scores = pd.DataFrame(index=clustering_results.columns, columns=categorical_metadata.columns)

# # For each clustering result
# for clustering in clustering_results.columns:
#     # For each categorical metadata
#     for metadata in categorical_metadata.columns:
#         # Calculate the scores
#         ami = metrics.adjusted_mutual_info_score(categorical_metadata[metadata], clustering_results[clustering])
#         ari = metrics.adjusted_rand_score(categorical_metadata[metadata], clustering_results[clustering])
#         fmi = metrics.fowlkes_mallows_score(categorical_metadata[metadata], clustering_results[clustering])
#         homogeneity, completeness, v_measure = metrics.homogeneity_completeness_v_measure(categorical_metadata[metadata], clustering_results[clustering])

#         # Store the scores in the corresponding DataFrame
#         ami_scores.loc[clustering, metadata] = ami
#         ari_scores.loc[clustering, metadata] = ari
#         fmi_scores.loc[clustering, metadata] = fmi
#         homogeneity_scores.loc[clustering, metadata] = homogeneity
#         completeness_scores.loc[clustering, metadata] = completeness
#         v_measure_scores.loc[clustering, metadata] = v_measure

# # Save the DataFrames as CSV files
# ami_scores.to_csv(ami_path)
# ari_scores.to_csv(ari_path)
# fmi_scores.to_csv(fmi_path)
# homogeneity_scores.to_csv(homogeneity_path)
# completeness_scores.to_csv(completeness_path)
# v_measure_scores.to_csv(v_path)