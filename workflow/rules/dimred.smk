####### perform Principal Component Analysis (PCA) #######
rule pca:
    input:
        unpack(get_sample_paths),
    output:
        result_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_object.pickle'),
        result_data = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_data.csv'),
        result_data_small = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_data_small.csv'),
        result_loadings = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_loadings.csv'),
        result_loadings_small = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_loadings_small.csv'),
        result_var = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_var.csv'),
        result_axes = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_{parameters}_axes.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","PCA_{sample}_{parameters}.log"),
    params:
        partition = config.get("partition"),
        samples_by_features = get_data_orientation,
    script:
        "../scripts/pca.py"
        
        
####### perform Uniform Manifold Approximation and Projection (UMAP) #######

# generate parametrized knn graphs using the UMAP package
rule umap_graph:
    input:
        unpack(get_sample_paths),
    output:
        result_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','UMAP_{metric}_{n_neighbors}_graph.pickle'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","umap_{sample}_{metric}_{n_neighbors}.log"),
    params:
        partition=config.get("partition"),
        samples_by_features = get_data_orientation,
        metric = lambda w: "{}".format(w.metric),
        n_neighbors = lambda w: "{}".format(w.n_neighbors),
    script:
        "../scripts/umap_graph.py"
        
        
# embed parametrized knn graphs using the UMAP package into lower dimensional space
rule umap_embed:
    input:
        get_umap_sample_paths,
    output:
        result_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','UMAP_{metric}_{n_neighbors}_{min_dist}_{n_components}_object.pickle'),
        result_data = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','UMAP_{metric}_{n_neighbors}_{min_dist}_{n_components}_data.csv'),
        result_axes = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','UMAP_{metric}_{n_neighbors}_{min_dist}_{n_components}_axes.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","umap_{sample}_{metric}_{n_neighbors}_{min_dist}_{n_components}.log"),
    params:
        partition=config.get("partition"),
        samples_by_features = get_data_orientation,
        metric = lambda w: "{}".format(w.metric),
        n_neighbors = lambda w: "{}".format(w.n_neighbors),
        min_dist = lambda w: "{}".format(w.min_dist),
        n_components = lambda w: "{}".format(w.n_components),
        densmap = 0,
    script:
        "../scripts/umap_embed.py"
        
# embed parametrized knn graphs using the UMAP package with densMAP flag into lower dimensional space
rule densmap_embed:
    input:
        get_umap_sample_paths,
    output:
        result_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','densMAP','densMAP_{metric}_{n_neighbors}_{min_dist}_{n_components}_object.pickle'),
        result_data = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','densMAP','densMAP_{metric}_{n_neighbors}_{min_dist}_{n_components}_data.csv'),
        result_axes = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','densMAP','densMAP_{metric}_{n_neighbors}_{min_dist}_{n_components}_axes.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","densmap_{sample}_{metric}_{n_neighbors}_{min_dist}_{n_components}.log"),
    params:
        partition=config.get("partition"),
        samples_by_features = get_data_orientation,
        metric = lambda w: "{}".format(w.metric),
        n_neighbors = lambda w: "{}".format(w.n_neighbors),
        min_dist = lambda w: "{}".format(w.min_dist),
        n_components = lambda w: "{}".format(w.n_components),
        densmap = 1,
    script:
        "../scripts/umap_embed.py"