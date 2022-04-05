# PCA scatter plot panel by metadata
rule plot_pca_metadata:
    input:
        unpack(get_pca_paths),
    output:
        metadata_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','plots','PCA_metadata.png'),
                               caption="../report/pca_2d_metadata.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "32000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_metadata_{sample}_PCA.log"),
    params:
        partition=config.get("partition"),
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"]
    script:
        "../scripts/plot_2d.R"
        
# dimred scatter plot panel by metadata
rule plot_dimred_metadata:
    input:
        unpack(get_dimred_paths),
    output:
        metadata_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','plots','{method}_{parameters}_metadata.png'),
                               caption="../report/dimred_2d_metadata.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem_small", "8000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_metadata_{sample}_{method}_{parameters}.log"),
    params:
        partition=config.get("partition"),
        size = config["scatterplot2d"]["size"],
        alpha = config["scatterplot2d"]["alpha"]
    script:
        "../scripts/plot_2d.R"

# PCA scree plot and cumulative variance
rule plot_pca_diagnostics:
    input:
        var_data = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','PCA_var.csv'),
    output:
        diagnostics_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','PCA','plots','PCA_diagnostics.png'),
                               caption="../report/pca_diagnostics.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem_small", "8000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","plot_diagnostics_{sample}_PCA.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_pca.R"
        
        
        
# plot UMAP diagnostic visualizations
rule plot_umap_diagnostics:
    input:
        umap_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','UMAP_{parameters}_object.pickle'),
    output:
        diagnostics_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','plots','UMAP_{parameters}_diagnostics.png'),
                               caption="../report/umap_diagnostics.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","plot_diagnostics_{sample}_UMAP_{parameters}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_umap_diagnostics.py"
        
        
# plot UMAP connectivity visualizations
rule plot_umap_connectivity:
    input:
        umap_object = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','UMAP_{parameters}_object.pickle'),
    output:
        connectivity_plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','UMAP','plots','UMAP_{parameters}_connectivity.png'),
                               caption="../report/umap_connectivity.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/umap.yaml"
    log:
        os.path.join("logs","rules","plot_connectivity_{sample}_UMAP_{parameters}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/plot_umap_connectivity.py"