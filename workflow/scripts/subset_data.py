#### subset data for usage as metadata in subsequent feature plots ####

#### libraries
# general
import os
import pandas as pd

#### configurations

# ipnuts
data_path = snakemake.input["data"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/condition_24h_cytokines/counts/CORRECTED_RNA.csv"
# outputs
result_data_path = snakemake.output[0] #"/nobackup/lab_bock/projects/macroIC/results/Lee2020NatGenet/unsupervised_analysis/merged_NORMALIZED/metadata_features.csv"

result_dir = os.path.dirname(result_data_path)

# parameters
samples_by_features = int(snakemake.params['samples_by_features']) # 0
features_to_plot = set(snakemake.params["features_to_plot"]) # set(['FCN1', 'TGFBR1', 'TGFBR2', 'TNFRSF1A','IL6R', 'IFNGR1', 'IFNGR2', 'IFNAR1', 'IFNG', 'IFNB1', 'IL6', 'TNF', 'TGFB1', 'TGFB2'])

# make directory if not existing
if not os.path.exists(result_dir):
    os.makedirs(result_dir, exist_ok=True)

### load data

# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T
    
### check overlap with columns & subset data
features_to_plot = list(features_to_plot.intersection(set(data.columns)))

if len(features_to_plot)!=0:
    data = data.loc[:,features_to_plot]
else:
    print("requested features to plot are not in the provided data, first 10 features will be plotted instead")
    data = data.iloc[:,:10]

# save data
data.to_csv(result_data_path)