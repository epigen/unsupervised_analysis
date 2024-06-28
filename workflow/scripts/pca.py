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
data_path = snakemake.input["data"]

# outputs
result_object_path = snakemake.output["result_object"]
result_data_path = snakemake.output["result_data"]
result_data_small_path = snakemake.output["result_data_small"]
result_loadings_path = snakemake.output["result_loadings"]
result_loadings_small_path = snakemake.output["result_loadings_small"]
result_var_path = snakemake.output["result_var"]
result_axes_path = snakemake.output["result_axes"]

# parameters
samples_by_features = int(snakemake.params['samples_by_features'])
n_components = snakemake.config["pca"]["n_components"]
svd_solver = snakemake.config["pca"]["svd_solver"]

### load data

# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T
    
### transform data
    
# unsupervised PCA
pca_obj = PCA(n_components=n_components,
              copy=True, 
              whiten=False, 
              svd_solver=svd_solver,
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
data_df.iloc[:,:min(10,data_df.shape[1])].to_csv(result_data_small_path)

# save loadings
loadings = pd.DataFrame(pca_obj.components_.T, columns = data_df.columns, index=data.columns)
loadings.to_csv(result_loadings_path)
loadings.iloc[:,:min(10,data_df.shape[1])].to_csv(result_loadings_small_path)

# save explained variance
axes_info_df = pd.DataFrame(pca_obj.explained_variance_ratio_)
axes_info_df.to_csv(result_var_path)

# save axes information for visualization
axes_info_df.columns = ['label']
axes_info_df['label'] = ["PC{} ({}%)".format(idx+1, round(100*var,2)) for idx, var in axes_info_df['label'].items()]
axes_info_df.to_csv(result_axes_path)
