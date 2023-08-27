# perform Leiden clustering
rule leiden_cluster:
    input:
        get_umap_sample_paths,
    output:
        clustering = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','Leiden','Leiden_{metric}_{n_neighbors}_{partition_type}_{resolution}_clustering.txt'),
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
        
# aggregate clustering results
rule aggregate_clustering_results:
    input:
        get_clustering_paths,
    output:
        aggregated_clusterings = os.path.join(config["result_path"],'unsupervised_analysis','{sample}','{method}','{method}_clusterings.csv'),
    log:
        os.path.join("logs","rules","aggregate_clustering_results_{sample}_{method}.log"),
    run:
        # dictionary to hold the data
        agg_clust = {}

        # read each clustering result and add to data dict
        for filename in input[1:]:
            clust_tmp = pd.read_csv(filename, header=None).squeeze("columns")
            agg_clust[os.path.splitext(os.path.basename(filename))[0]] = clust_tmp

        # Load the index from the metadata
        metadata_df = pd.read_csv(input[0], index_col=0)
        
        # convert the dictionary to a DataFrame
        agg_clust_df = pd.DataFrame(agg_clust, index=metadata_df.index)
        
        # Write the DataFrame to a CSV file
        agg_clust_df.to_csv(output.aggregated_clusterings, index=True)

