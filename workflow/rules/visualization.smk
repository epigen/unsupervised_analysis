########## FEATURE PLOTS ##########

# prepare "metadata" for feature plot
rule prep_feature_plot:
    input:
        unpack(get_sample_paths),
    output:
        metadata_features = os.path.join(result_path,'{sample}','metadata_features.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap_leiden.yaml"
    log:
        os.path.join("logs","rules","prep_feature_plot_{sample}.log"),
    params:
        partition = config.get("partition"),
        samples_by_features = get_data_orientation,
        features_to_plot = config["features_to_plot"],
    script:
        "../scripts/subset_data.py"
        
# dimred feature scatter plots
rule plot_dimred_features:
    input:
        unpack(get_dimred_features_paths),
    output:
        plot = report(directory(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}_{n_components}','features')),
                      patterns=["{feature}.png"],
                      caption="../report/dimred_2d_features.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{sample}",
                      labels={
                          "method": "{method}",
                          "parameters": "{parameters}",
                          "dimensions": "{n_components}",
                          "type": "features",
                          "content": "{feature}",
                      }
                     ),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_features_{sample}_{method}_{parameters}_{n_components}.log"),
    params:
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"],
        partition=config.get("partition"),
    script:
        "../scripts/plot_2d.R"

########## METADATA PLOTS ##########
        
# dimred scatter plot panel by metadata
rule plot_dimred_metadata:
    input:
        unpack(get_dimred_paths),
    output:
        plot = report(directory(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}_{n_components}','metadata')),
                      patterns=["{metadata}.png"],
                      caption="../report/dimred_2d_metadata.rst",
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{sample}",
                      labels={
                          "method": "{method}",
                          "parameters": "{parameters}",
                          "dimensions": "{n_components}",
                          "type": "metadata",
                          "content": "{metadata}",
                      }
                     ),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_metadata_{sample}_{method}_{parameters}_{n_components}.log"),
    params:
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"],
        partition=config.get("partition"),
    script:
        "../scripts/plot_2d.R"

########## DIAGNOSTIC PLOTS ##########
        
# PCA scree plot, cumulative variance plot, pairs plot, and loadings plot
rule plot_pca_diagnostics:
    input:
        unpack(get_dimred_paths)
    output:
        diagnostics_plot = report(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}','diagnostics','variance.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_{}".format(config["project_name"], module_name),
                               subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{parameters}",
                                  "dimensions": "-",
                                  "type": "diagnostics",
                                  "content": "variance",
                      }),
        pairs_plot = report(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}','diagnostics','pairs.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_{}".format(config["project_name"], module_name),
                               subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{parameters}",
                                  "dimensions": "-",
                                  "type": "diagnostics",
                                  "content": "pairs",
                      }),
        loadings_plot = report(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}','diagnostics','loadings.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_{}".format(config["project_name"], module_name),
                               subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{parameters}",
                                  "dimensions": "-",
                                  "type": "diagnostics",
                                  "content": "loadings",
                      }),
        loadings_lollipop_plot = report(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}','diagnostics','loadings_lollipop.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_{}".format(config["project_name"], module_name),
                               subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{parameters}",
                                  "dimensions": "-",
                                  "type": "diagnostics",
                                  "content": "loadings lollipop",
                      }),
    resources:
        mem_mb=config.get("mem", "8000"),
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
        umap_object = os.path.join(result_path,'{sample}','{method}','{method}_{parameters}_object.pickle'),
    output:
        diagnostics_plot = report(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}','diagnostics','diagnostics.png'),
                               caption="../report/umap_diagnostics.rst", 
                               category="{}_{}".format(config["project_name"], module_name),
                               subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{parameters}",
                                  "dimensions": "-",
                                  "type": "diagnostics",
                                  "content": "diagnostics",
                      }),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap_leiden.yaml"
    log:
        os.path.join("logs","rules","plot_diagnostics_{sample}_{method}_{parameters}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_umap_diagnostics.py"
        
        
# plot UMAP & densMAP connectivity visualizations
rule plot_umap_connectivity:
    input:
        umap_object = os.path.join(result_path,'{sample}','{method}','{method}_{parameters}_object.pickle'),
    output:
        connectivity_plot = report(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}','diagnostics','connectivity.png'),
                               caption="../report/umap_connectivity.rst", 
                               category="{}_{}".format(config["project_name"], module_name),
                               subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{parameters}",
                                  "dimensions": "-",
                                  "type": "diagnostics",
                                  "content": "connectivity",
                      }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap_leiden.yaml"
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
        plot = os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}_{n_components}','interactive.html'),
    resources:
        mem_mb=config.get("mem", "8000"),
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
        unpack(get_heatmap_paths),
    output:
        plot = report(os.path.join(result_path,'{sample}','Heatmap','plots','Heatmap_{metric}_{method}.png'),
                      caption="../report/heatmap.rst",
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{sample}",
                                 labels={
                                  "method": "{method}",
                                  "parameters": "{metric}",
                                  "dimensions": "-",
                                  "type": "Heatmap",
                                  "content": "Heatmap",
                      }),
    resources:
        # dynamic memory allocation based on input size and attempts (multiple attempts can be triggered with --retries X)
        mem_mb=lambda wildcards, attempt: attempt*int(config.get("mem", "16000")),#lambda wildcards, input, attempt: max(int(config.get("mem", "16000")),((input.size//1000000) * attempt * 70)),#config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ComplexHeatmap.yaml"
    log:
        os.path.join("logs","rules","plot_heatmap_{sample}_{metric}_{method}.log"),
    params:
        partition = config.get("partition"),
        samples_by_features = get_data_orientation,
    script:
        "../scripts/plot_heatmap.R"
        
        
########## CLUSTERING PLOTS ##########
        
# dimred 2D scatter plot panel by clustering
rule plot_dimred_clustering:
    input:
        unpack(get_metadata_clustering_paths),
    output:
        plot = report(directory(os.path.join(result_path,'{sample}','{method}','plots','{method}_{parameters}_{n_components}','clustering')),
                      patterns=["{clustering}.png"],
                      caption="../report/dimred_2d_clusterings.rst", 
                      category="{}_{}".format(config["project_name"], module_name),
                      subcategory="{sample}",
                      labels={
                          "method": "{method}",
                          "parameters": "{parameters}",
                          "dimensions": "{n_components}",
                          "type": "clustering",
                          "content": "{clustering}",
                      }
                     ),
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
