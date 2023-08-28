########## FEATURE PLOTS ##########

# prepare "metadata" for feature plot
rule prep_feature_plot:
    input:
        unpack(get_sample_paths),
    output:
        os.path.join(config["result_path"],'unsupervised_analysis','{sample}','metadata_features.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","prep_feature_plot_{sample}.log"),
    params:
        partition = config.get("partition"),
        samples_by_features = get_data_orientation,
        features_to_plot = config["features_to_plot"],
    script:
        "../scripts/subset_data.py"
        
# dimred scatter plot panel by features
rule plot_dimred_features:
    input:
        unpack(get_dimred_features_paths),
    output:
        plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_{n_components}_features.png'),
                               caption="../report/dimred_2d_features.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_features_{sample}_{method}_{parameters}_{n_components}.log"),
    params:
        partition=config.get("partition"),
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"]
    script:
        "../scripts/plot_2d.R"

########## METADATA PLOTS ##########
        
# dimred scatter plot panel by metadata
rule plot_dimred_metadata:
    input:
        unpack(get_dimred_paths),
    output:
        plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_{n_components}_metadata.png'),
                               caption="../report/dimred_2d_metadata.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_metadata_{sample}_{method}_{parameters}_{n_components}.log"),
    params:
        partition=config.get("partition"),
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"]
    script:
        "../scripts/plot_2d.R"

########## DIAGNOSTIC PLOTS ##########
        
# PCA scree plot, cumulative variance plot, pairs plot, and loadings plot
rule plot_pca_diagnostics:
    input:
        unpack(get_dimred_paths)
    output:
        diagnostics_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_variance.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
        pairs_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_pairs.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
        loadings_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_loadings.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem_small", "8000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_{method}_diagnostics_{sample}_{parameters}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_pca.R"
        

# plot UMAP & densMAP diagnostic visualizations
rule plot_umap_diagnostics:
    input:
        umap_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','{method}_{parameters}_object.pickle'),
    output:
        diagnostics_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_diagnostics.png'),
                               caption="../report/umap_diagnostics.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","plot_diagnostics_{sample}_{method}_{parameters}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_umap_diagnostics.py"
        
        
# plot UMAP & densMAP connectivity visualizations
rule plot_umap_connectivity:
    input:
        umap_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','{method}_{parameters}_object.pickle'),
    output:
        connectivity_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_connectivity.png'),
                               caption="../report/umap_connectivity.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","plot_connectivity_{sample}_{method}_{parameters}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_umap_connectivity.py"

########## INTERACTIVE PLOTS ##########

# plot interactive HTML plots of 2D and 3D embeddings using plotly
rule plot_dimred_interactive:
    input:
        unpack(get_dimred_paths),
    output:
        plot = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_{n_components}_interactive.html'),
    resources:
        mem_mb=config.get("mem_small", "8000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/plotly.yaml"
    log:
        os.path.join("logs","rules","plot_interactive_{sample}_{method}_{parameters}_{n_components}.log"),
    params:
        partition=config.get("partition"),
        n_components = lambda w: "{}".format(w.n_components),
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"]
    script:
        "../scripts/plot_interactive.py"
        
        
########## HEATMAP PLOTS ##########

rule plot_heatmap:
    input:
        unpack(get_sample_paths),
    output:
        plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','Heatmap','plots','Heatmap_{method}_{metric}.png'),
                      caption="../report/heatmap.rst",
                      category="{}_unsupervised_analysis".format(config["project_name"]),
                      subcategory="{sample}"),
    resources:
        # dynamic memory allocation based on input size and attempts (multiple attempts can be triggered with --retries X)
        mem_mb=lambda wildcards, attempt: attempt*int(config.get("mem", "16000")),#lambda wildcards, input, attempt: max(int(config.get("mem", "16000")),((input.size//1000000) * attempt * 70)),#config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ComplexHeatmap.yaml"
    log:
        os.path.join("logs","rules","plot_heatmap_{sample}_{method}_{metric}.log"),
    params:
        partition = config.get("partition"),
        samples_by_features = get_data_orientation,
        metric = lambda w: "{}".format(w.metric),
        cluster_method = lambda w: "{}".format(w.method)
    script:
        "../scripts/plot_heatmap.R"
        
        
########## CLUSTERING PLOTS ##########
        
# dimred 2D scatter plot panel by clustering
rule plot_dimred_clustering:
    input:
        unpack(get_metadata_clustering_paths),
    output:
        plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_{n_components}_clustering.png'),
                               caption="../report/dimred_2d_clusterings.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_clustering_{sample}_{method}_{parameters}_{n_components}.log"),
    params:
        partition=config.get("partition"),
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"]
    script:
        "../scripts/plot_2d.R"