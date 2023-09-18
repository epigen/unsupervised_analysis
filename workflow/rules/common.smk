##### utility functions #####

########## GENERAL ##########

def get_sample_paths(wildcards):
    return {'data': annot.loc[wildcards.sample,'data'],
           'metadata': annot.loc[wildcards.sample,"metadata"]
           }

def get_data_orientation(wildcards):
    return int(annot.loc[wildcards.sample,'samples_by_features'])

########## DIMENSIONALITY REDUCTION ##########

def get_umap_sample_paths(wildcards):
    return [annot.loc[wildcards.sample,'data'],
           os.path.join(config["result_path"],'unsupervised_analysis','{}'.format(wildcards.sample),'UMAP','UMAP_{}_'.format(wildcards.metric)+'{}'.format(max(config["umap"]["n_neighbors"]))+'_graph.pickle')]

def get_dimred_paths(wildcards):
    path_dict = {}
    
    if wildcards.method=="PCA":
        path_dict['dimred_data'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_data_small.csv'.format(wildcards=wildcards))
        path_dict['dimred_axes'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_axes.csv'.format(wildcards=wildcards))
        path_dict['dimred_var'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_var.csv'.format(wildcards=wildcards))
        path_dict['dimred_loadings'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_loadings_small.csv'.format(wildcards=wildcards))
#         return {
#             'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_data_small.csv'.format(wildcards=wildcards)),
#             'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_axes.csv'.format(wildcards=wildcards)),
#             'dimred_var': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_var.csv'.format(wildcards=wildcards)),
#             'dimred_loadings': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_loadings_small.csv'.format(wildcards=wildcards)),
#             'metadata': annot.loc[wildcards.sample,"metadata"],
#             'metadata_features': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
#                }
    else:
        path_dict['dimred_data'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_data.csv'.format(wildcards=wildcards))
        path_dict['dimred_axes'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_axes.csv'.format(wildcards=wildcards))
#         return {
#             'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_data.csv'.format(wildcards=wildcards)),
#             'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_axes.csv'.format(wildcards=wildcards)),
#             'metadata': annot.loc[wildcards.sample,"metadata"],
#             'metadata_features': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
#                }
    
    # add metadata
    path_dict['metadata'] = annot.loc[wildcards.sample,"metadata"]
    # add features
    path_dict['metadata_features'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
    # add clustering results
    if len(cluster_methods) > 0:
        path_dict['metadata_clusterings'] = os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_clusterings.csv')
    
    return path_dict
    
def get_dimred_features_paths(wildcards):
    
    if wildcards.method=="PCA":
        return {
            'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_data_small.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_axes.csv'.format(wildcards=wildcards)),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
               }
    else:
        return {
            'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_axes.csv'.format(wildcards=wildcards)),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
               }


########## CLUSTERING ##########

# get paths for clustification
def get_clustification_paths(wildcards):
    return [annot.loc[wildcards.sample,'data'],
            os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'{}'.format(config["clustification"]["method"]),'{}_clusterings.csv'.format(config["clustification"]["method"]))
           ]

# get all clustering results of one method to be aggregated into {method}/{method}_clusterings.csv
def get_clustering_paths(wildcards):
    path_list = [annot.loc[wildcards.sample,"metadata"]]
    
    if wildcards.method=="Leiden":
        leiden_parameters = []
        
        # add resolution parameter only to relevant partition_type algorithms, otherwise NA
        for partition_type in config["leiden"]["partition_types"]:
            if partition_type in ["RBConfigurationVertexPartition", "RBERVertexPartition", "CPMVertexPartition"]:
                leiden_parameters = leiden_parameters + ["{}_{}".format(partition_type, res) for res in config["leiden"]["resolutions"]]
            else:
                leiden_parameters.append("{}_NA".format(partition_type))
        
        path_list = path_list + expand(os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'Leiden','Leiden_{metric}_{n_neighbors}_{leiden_parameters}_clustering.csv'),
                                metric=config["leiden"]["metrics"],
                                n_neighbors=config["leiden"]["n_neighbors"],
                                leiden_parameters=leiden_parameters,
#                                 partition_type=config["leiden"]["partition_types"],
#                                 resolution=config["leiden"]["resolutions"]
                                   )
    return path_list

# get all aggregated clustering results across methods to be aggregated into {sample}/metadata_clusterings.csv
def get_aggregated_clustering_paths(wildcards):
    return expand(os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'{method}','{method}_clusterings.csv'), method=cluster_methods)
    
# get the aggregated clustering results acros methods for visualization
def get_metadata_clustering_paths(wildcards):
    return {
            'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_axes.csv'.format(wildcards=wildcards)),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample, "metadata_clusterings.csv")
    }

########## CLUSTER VALIDATION ##########

# get the input paths for clustree analysis depending on requested content type
def get_clustree_paths(wildcards):
    
    if wildcards.content=="features":
        return {
            'metadata_clustering': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample, "metadata_clusterings.csv"),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
        }
    else:
        return {
            'metadata_clustering': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample, "metadata_clusterings.csv"),
            'metadata': annot.loc[wildcards.sample,"metadata"]
        }

# get paths to determine external cluster indices
def get_external_validation_paths(wildcards):
    return {'clusterings': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample, "metadata_clusterings.csv"),
            'metadata': annot.loc[wildcards.sample,"metadata"]
           }

# for plotting heatmaps of cluster indices
def get_validation_paths(wildcards):
    return {
        idx: os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample, "cluster_validation", "external_index_{}.csv".format(idx)) for idx in indices_external
    }



