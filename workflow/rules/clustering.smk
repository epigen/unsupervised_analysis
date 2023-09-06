# perform Leiden clustering
rule leiden_cluster:
    input:
        get_umap_sample_paths,
    output:
        clustering = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','Leiden','Leiden_{metric}_{n_neighbors}_{partition_type}_{resolution}_clustering.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/leiden.yaml"
    log:
        os.path.join("logs","rules","leiden_{sample}_{metric}_{n_neighbors}_{partition_type}_{resolution}_clustering.log"),
    params:
        partition=config.get("partition"),
        samples_by_features = get_data_orientation,
        metric = lambda w: "{}".format(w.metric),
        n_neighbors = lambda w: "{}".format(w.n_neighbors),
        partition_type = lambda w: "{}".format(w.partition_type),
        resolution = lambda w: "{}".format(w.resolution),
        n_iterations = config["leiden"]["n_iterations"]
    script:
        "../scripts/leiden_cluster.py"
        
# perform clustification based on initial clustering
rule clustification:
    input:
        get_clustification_paths,
    output:
        clustering = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','clustification','clustification_clusterings.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","clustification_{sample}_clusterings.log"),
    params:
        partition=config.get("partition"),
        samples_by_features = get_data_orientation,
    script:
        "../scripts/clustification.py"
        
# aggregate clustering results per method
rule aggregate_clustering_results:
    input:
        get_clustering_paths,
    output:
        aggregated_clusterings = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','{method}_clusterings.csv'),
    log:
        os.path.join("logs","rules","aggregate_clustering_results_{sample}_{method}.log"),
    run:
        # list to hold the individual clusterings
        agg_clust = []

        # read each clustering result and add to list
        for filename in input[1:]:
            clust_tmp = pd.read_csv(filename, header=0, index_col=0)#.squeeze("columns")
            agg_clust.append(clust_tmp)

        # convert list to dataframe
        agg_clust_df = pd.concat(agg_clust, axis=1)
        
        # Write the DataFrame to a CSV file
        agg_clust_df.to_csv(output.aggregated_clusterings, index=True)
        
# aggregate clustering results across methods
rule aggregate_all_clustering_results:
    input:
        get_aggregated_clustering_paths,
    output:
        metadata_clusterings = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','metadata_clusterings.csv'),
    log:
        os.path.join("logs","rules","aggregate_all_clustering_results_{sample}.log"),
    run:
        # list to hold the data
        agg_clust = []

        # read each clustering result and add to data dict
        for filename in input:
            agg_clust.append(pd.read_csv(filename, header=0, index_col=0))

        
        # convert the dictionary to a DataFrame
        agg_clust_df = pd.concat(agg_clust, axis=1)
        
        # Write the DataFrame to a CSV file
        agg_clust_df.to_csv(output.metadata_clusterings, index=True)
        
# clustree analysis
rule clustree_analysis:
    input:
        unpack(get_clustree_paths),
    output:
        plot = report(os.path.join(config["result_path"],'unsupervised_analysis','{sample}','clustree','clustree_{content}.png'),
                               caption="../report/clustree.rst", 
                               category="{}_unsupervised_analysis".format(config["project_name"]), 
                               subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/clustree.yaml"
    log:
        os.path.join("logs","rules","clustree_{sample}_{content}.log"),
    params:
        partition=config.get("partition"),
        content = lambda w: "{}".format(w.content),
        count_filter = config["clustree"]["count_filter"],
        prop_filter = config["clustree"]["prop_filter"],
        layout = config["clustree"]["layout"],
        categorical_label_option = config["clustree"]["categorical_label_option"],
        numerical_aggregation_option = config["clustree"]["numerical_aggregation_option"],
        custom_metadata = config["clustree"]["custom_metadata"],
    script:
        "../scripts/clustree.R"