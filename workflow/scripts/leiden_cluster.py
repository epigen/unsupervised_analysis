#### load knn graph from UMAP (pyNNDescent) and perform Leiden clustering on the graph ####

#### libraries
import os
import pickle
import numpy as np
# Leiden algorithm
import leidenalg as la
from igraph import Graph
from scipy.sparse import csr_matrix
import pandas as pd

# helper function adapted from here: https://dynamo-release.readthedocs.io/en/latest/_modules/dynamo/tools/connectivity.html
def knn_to_adj(knn_indices: np.ndarray, knn_weights: np.ndarray) -> csr_matrix:
    """Convert a knn graph's indices and weights to an adjacency matrix of the corresponding nearest neighbor graph.

    Args:
        knn_indices: the matrix (n x n_neighbors) storing the indices for each node's n_neighbors nearest neighbors in
            the knn graph.
        knn_weights: the matrix (n x n_neighbors) storing the weights on the edges for each node's n_neighbors nearest
            neighbors in the knn graph.

    Returns:
        The converted adjacency matrix (n x n) of the corresponding nearest neighbor graph.
    """

    adj = csr_matrix(
        (
            knn_weights.flatten(),
            (
                np.repeat(knn_indices[:, 0], knn_indices.shape[1]),
                knn_indices.flatten(),
            ),
        )
    )
    adj.eliminate_zeros()

    return adj

#### configurations

# inputs
data_path = os.path.join(snakemake.input[0]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/data/digits_data.csv"
graph_path = os.path.join(snakemake.input[1]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/UMAP/UMAP_euclidean_100_graph.pickle"

# outputs
result_path = os.path.join(snakemake.output["clustering"]) # "/research/home/sreichl/projects/unsupervised_analysis/.test/results/unsupervised_analysis/digits/Leiden/Leiden_RBConfigurationVertexPartition_1.csv"

# UMAP parameters for small data (<11 observations)
samples_by_features = int(snakemake.params['samples_by_features']) #0
metric = snakemake.params['metric'] #"correlation"
n_neighbors = int(snakemake.params['n_neighbors']) #100

# Leiden algorithm parameters
n_iterations = int(snakemake.params["n_iterations"]) # 2
# Get the partition method from the leidenalg module
partition_type = getattr(la, str(snakemake.params["partition_type"])) #"RBConfigurationVertexPartition"
# Get the kwargs from the config dict
if str(snakemake.params["partition_type"]) in ["RBConfigurationVertexPartition", "RBERVertexPartition", "CPMVertexPartition"]:
    la_kwargs = {"resolution_parameter": float(snakemake.params["resolution"])} #0.05 #1.0
else:
    la_kwargs = {}

result_dir = os.path.dirname(result_path)
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
    # run UMAP
    umap_obj = umap.umap_.UMAP(n_neighbors=n_neighbors, 
                         n_components=2, 
                         metric=metric, 
                         metric_kwds=None, 
                         output_metric='euclidean', 
                         output_metric_kwds=None, 
                         n_epochs=None, 
                         learning_rate=1.0, 
                         init='spectral', 
                         min_dist=0.1, 
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
                         precomputed_knn=(None, None, None)
                        ).fit(data)
    # extract graph
    adj_coo = umap_obj.graph_.tocoo()
else:
    # load pre-computed KNN graph
    with open(graph_path, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        knn = pickle.load(f)
    # convert to CSR (Compressed Sparse Row) format
    adj_csr = knn_to_adj(knn[0], knn[1])
    # Convert the CSR matrix to COO format
    adj_coo = adj_csr.tocoo()

# Convert the adjacency matrix in COO format to igraph object (weighted undirected graph)
edges = np.column_stack((adj_coo.row, adj_coo.col))
graph = Graph(edges.tolist(), directed=False)
graph.es['weight'] = adj_coo.data.tolist()


# Perform Leiden clustering
partition = la.find_partition(graph=graph, 
                              partition_type=partition_type,
                              initial_membership=None, # default
                              weights='weight',  # default: None
                              n_iterations=n_iterations,
                              max_comm_size=0, 
                              seed=42,
                              **la_kwargs
                             )

# save clustering as CSV
clustering_name = os.path.splitext(os.path.basename(result_path))[0]
pd.DataFrame({clustering_name: partition.membership}, index=data.index).to_csv(result_path, index=True)

# with open(result_path, 'w') as f:
#     for label in partition.membership:
#         f.write(str(label) + '\n')
