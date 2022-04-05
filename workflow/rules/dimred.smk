# perform Principal Component Analysis (PCA)
rule pca:
    input:
        get_sample_paths,
    output:
        result_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_object.pickle'),
        result_data = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_data.csv'),
        result_var = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_var.csv'),
        result_axes = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_axes.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","pca_{sample}.log"),
    params:
        partition=config.get("partition"),
        samples_by_features = config["samples_by_features"],
    script:
        "../scripts/pca.py"
        
        
# perform Uniform Manifold Approximation and Projection (UMAP)

# generate parametrized knn graphs using the UMAP package
rule umap_graph:
    input:
        get_sample_paths,
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
        samples_by_features = config["samples_by_features"],
        metric = lambda w: "{}".format(w.metric),
        n_neighbors = lambda w: "{}".format(w.n_neighbors),
    script:
        "../scripts/umap_graph.py"
        
        
# embed parametrized knn graphs using the UMAP package in to lower dimensional space
rule umap_embed:
    input:
        get_umap_sample_paths,
#         graph_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP_{metric}_'+'{}'.format(max(config["umap"]["n_neighbors"]))+'_graph.pickle'),
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
        samples_by_features = config["samples_by_features"],
        metric = lambda w: "{}".format(w.metric),
        n_neighbors = lambda w: "{}".format(w.n_neighbors),
        min_dist = lambda w: "{}".format(w.min_dist),
        n_components = lambda w: "{}".format(w.n_components),
    script:
        "../scripts/umap_embed.py"