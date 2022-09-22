#### perform low dimensional embedding of the knn-graph using Uniform Manifold Approximation and Projection (UMAP) ####

#### libraries
# general
import os
import pickle
import pandas as pd
# dimensionality reduction
import umap

#### configurations

# ipnuts
data_path = snakemake.input[0] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/KOcall_NonTargeting/counts/CORRECTED_RNA.csv"
graph_object_path = snakemake.input[1] #snakemake.input["knn_object"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP_correlation_100_graph.pickle"
# outputs
result_object_path = snakemake.output["result_object"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP_correlation_100_0.1_2_object.pickle"
result_data_path = snakemake.output["result_data"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP_correlation_100_0.1_2_data.csv"
result_axes_path = snakemake.output["result_axes"] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/unsupervised_analysis/AKsmall_KOcall_NonTargeting_CORRECTED/UMAP_correlation_100_0.1_2_axes.csv"

result_dir = os.path.dirname(result_object_path)

# parameters
samples_by_features = int(snakemake.params['samples_by_features']) #0
metric = snakemake.params['metric'] #"correlation"
n_neighbors = int(snakemake.params['n_neighbors']) #100
min_dist = float(snakemake.params['min_dist']) #0.1
n_components = int(snakemake.params['n_components']) #2

# make directory if not existing
if not os.path.exists(result_dir):
    os.makedirs(result_dir, exist_ok=True)

### load data

# check data orientation to fit: samples/observations x features
if samples_by_features == 1:
    data = pd.read_csv(data_path, index_col=0)
else:
    data = pd.read_csv(data_path, index_col=0).T

# if less than 11 datapoints there is no pre-computed KNN graph
if data.shape[0]<11:
    print("no pre-computed KNN graph will be used")
    knn = (None, None, None)
else:
    # load pre-computed KNN graph
    with open(graph_object_path, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        knn = pickle.load(f)

### embed data in low dimensions
umap_obj=umap.umap_.UMAP(n_neighbors=n_neighbors, 
                         n_components=n_components, 
                         metric=metric, 
                         metric_kwds=None, 
                         output_metric='euclidean', 
                         output_metric_kwds=None, 
                         n_epochs=None, 
                         learning_rate=1.0, 
                         init='spectral', 
                         min_dist=min_dist, 
                         spread=1.0, 
                         low_memory=True, 
                         n_jobs=-1, 
                         set_op_mix_ratio=1.0, 
                         local_connectivity=1.0, 
                         repulsion_strength=1.0, 
                         negative_sample_rate=5, 
                         transform_queue_size=4.0, 
                         a=None, 
                         b=None, 
                         random_state=42, 
                         angular_rp_forest=False, 
                         target_n_neighbors=-1, 
                         target_metric='categorical', 
                         target_metric_kwds=None, 
                         target_weight=0.5, 
                         transform_seed=42, 
                         transform_mode='embedding', 
                         force_approximation_algorithm=False, 
                         verbose=False, 
                         tqdm_kwds=None, 
                         unique=False, 
                         densmap=False, 
                         dens_lambda=2.0, 
                         dens_frac=0.3, 
                         dens_var_shift=0.1, 
                         output_dens=False, 
                         disconnection_distance=None, 
                         precomputed_knn=knn
                        ).fit(data)

# data_embedding = umap_obj.fit_transform(data)

data_df = pd.DataFrame(umap_obj.embedding_, index=data.index,)
data_df = data_df.rename_axis(("sample_name"))
data_df.columns = ["UMAP_{}".format(str(idx+1)) for idx in data_df.columns]

### save data

# save umap object
with open(result_object_path, 'wb') as f:
    pickle.dump(umap_obj, f, pickle.HIGHEST_PROTOCOL)
    
# save transformed data
data_df.to_csv(result_data_path)

# save axes information for visualization
axes_info_df = pd.DataFrame(data_df.columns)
axes_info_df.columns = ['label']
axes_info_df['label'] = [label.replace("_", "")for label in axes_info_df['label']]
axes_info_df.to_csv(result_axes_path)

