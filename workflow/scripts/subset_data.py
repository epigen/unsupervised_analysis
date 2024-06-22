#### subset data for usage as metadata in subsequent feature plots ####

#### libraries
# general
import os
import pandas as pd

#### configurations

# ipnuts
data_path = snakemake.input["data"]
# outputs
metadata_features_path = snakemake.output["metadata_features"]

# parameters
samples_by_features = int(snakemake.params['samples_by_features'])
features_to_plot = set(snakemake.params["features_to_plot"])

### load data

# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T
    
### check if "ALL" features should be plotted and overlap with columns & subset data
if features_to_plot == {"ALL"}:
    features_to_plot = list(data.columns)
else:
    features_to_plot = list(features_to_plot.intersection(set(data.columns)))

# subset data
if len(features_to_plot)!=0:
    data = data.loc[:,features_to_plot]
else:
    print("requested features to plot are not in the provided data, first 10 features will be plotted instead")
    data = data.iloc[:,:10]

# save data
data.to_csv(metadata_features_path)