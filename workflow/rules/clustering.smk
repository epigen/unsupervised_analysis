# # perform Leiden clustering TODO
# rule clustering_leiden:
#     input:
#         ...,
#     output:
#         ...,
#     resources:
#         mem_mb=config.get("mem", "16000"),
#     threads: config.get("threads", 1)
#     conda:
#         "../envs/sklearn.yaml"
#     log:
#         os.path.join("logs","rules","pca_{sample}.log"),
#     params:
#         partition=config.get("partition"),
#         samples_by_features = config["samples_by_features"],
#     script:
#         "../scripts/pca.py"