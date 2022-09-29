#### perform Principal Component Analysis (PCA) using sklearn ####

#### libraries
# general
import os
import pickle
import pandas as pd
# dimensionality reduction
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#### configurations

# ipnuts
data_path = snakemake.input[0] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/condition_24h_cytokines/counts/CORRECTED_RNA.csv"
# outputs
result_object_path = snakemake.output["result_object"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_condition_24h_cytokines_CORRECTED/PCA_object.pickle"
result_data_path = snakemake.output["result_data"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_condition_24h_cytokines_CORRECTED/PCA_data.csv"
result_loadings_path = snakemake.output["result_loadings"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_condition_24h_cytokines_CORRECTED/PCA_loadings.csv"
result_var_path = snakemake.output["result_var"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_condition_24h_cytokines_CORRECTED/PCA_var.csv"
result_axes_path = snakemake.output["result_axes"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_condition_24h_cytokines_CORRECTED/PCA_axes.csv"

result_dir = os.path.dirname(result_object_path)

# parameters
samples_by_features = int(snakemake.params['samples_by_features']) #0

# make directory if not existing
if not os.path.exists(result_dir):
    os.makedirs(result_dir, exist_ok=True)

### load data

# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T
    
### transform data
    
# unsupervised PCA to yield all components
pca_obj = PCA(n_components=None, # all components are kept
              copy=True, 
              whiten=False, 
              svd_solver='full', # run exact full SVD
              tol=0.0, 
              iterated_power='auto', 
              random_state=42
             )

data_pca = pca_obj.fit_transform(StandardScaler().fit_transform(data))

data_df = pd.DataFrame(data_pca, index=data.index,)
data_df = data_df.rename_axis(("sample_name"))
data_df.columns = ["PC_{}".format(str(idx+1)) for idx in data_df.columns]

### save data

# save pca object
with open(result_object_path, 'wb') as f:
    pickle.dump(pca_obj, f, pickle.HIGHEST_PROTOCOL)
    
# save transformed data
data_df.to_csv(result_data_path)

# save loadings
loadings = pd.DataFrame(pca_obj.components_.T, columns = data_df.columns, index=data.columns)
loadings.to_csv(result_loadings_path)

# save explained variance
axes_info_df = pd.DataFrame(pca_obj.explained_variance_ratio_)
axes_info_df.to_csv(result_var_path)

# save axes information for visualization
axes_info_df.columns = ['label']
axes_info_df['label'] = ["PC{} ({}%)".format(idx+1, round(100*var,2)) for idx, var in axes_info_df['label'].iteritems()]
axes_info_df.to_csv(result_axes_path)
