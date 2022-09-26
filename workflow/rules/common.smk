##### utility functions #####

def get_sample_paths(wildcards):
    return [annot.loc[wildcards.sample,'data']]

def get_umap_sample_paths(wildcards):
    return [annot.loc[wildcards.sample,'data'],
           os.path.join(config["result_path"],'unsupervised_analysis','{}'.format(wildcards.sample),'UMAP','UMAP_{}_'.format(wildcards.metric)+'{}'.format(max(config["umap"]["n_neighbors"]))+'_graph.pickle')]

def get_pca_paths(wildcards):
    return {'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_axes.csv'.format(wildcards=wildcards)),
            'metadata': annot.loc[wildcards.sample,"metadata"]
           }

def get_pca_feature_paths(wildcards):
    return {'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'PCA','PCA_axes.csv'.format(wildcards=wildcards)),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv'.format(wildcards=wildcards))
           }

def get_dimred_paths(wildcards):
    return {'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_axes.csv'.format(wildcards=wildcards)),
            'metadata': annot.loc[wildcards.sample,"metadata"]
       }

def get_dimred_feature_paths(wildcards):
    return {'dimred_data': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_data.csv'.format(wildcards=wildcards)),
            'dimred_axes': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,wildcards.method,'{wildcards.method}_{wildcards.parameters}_axes.csv'.format(wildcards=wildcards)),
            'metadata': os.path.join(config["result_path"],'unsupervised_analysis',wildcards.sample,'metadata_features.csv'.format(wildcards=wildcards))
       }

def get_data_orientation(wildcards):
    return int(annot.loc[wildcards.sample,'samples_by_features'])