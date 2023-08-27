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
    
    if wildcards.method=="PCA":
        return {
            'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_data_small.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_axes.csv'.format(wildcards=wildcards)),
            'dimred_var': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_var.csv'.format(wildcards=wildcards)),
            'dimred_loadings': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_{wildcards.parameters}_loadings_small.csv'.format(wildcards=wildcards)),
            'metadata': annot.loc[wildcards.sample,"metadata"],
            'metadata_features': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
               }
    else:
        return {
            'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_axes.csv'.format(wildcards=wildcards)),
            'metadata': annot.loc[wildcards.sample,"metadata"],
            'metadata_features': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv')
               }
    
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

def get_clustering_paths(wildcards):
    if wildcards.method=="Leiden":
        return [annot.loc[wildcards.sample,"metadata"]] + expand(os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'Leiden','Leiden_{metric}_{n_neighbors}_{partition_type}_{resolution}_clustering.txt'),
                                metric=config["leiden"]["metrics"],
                                n_neighbors=config["leiden"]["n_neighbors"],
                                partition_type=config["leiden"]["partition_types"],
                                resolution=config["leiden"]["resolutions"]
                                   )
    else:
        return {}
    
# once there are multiple clustering methods an intermediate rule has to aggregate the results and save e.g., {sample}/all_clusterings.csv -> input for plotting
def get_aggregated_clustering_paths(wildcards):
    return {
            'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_{wildcards.n_components}_axes.csv'.format(wildcards=wildcards)),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,"Leiden", "Leiden_clusterings.csv")
    }