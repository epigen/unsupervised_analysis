# conda environment export rule to document the exact versions and builds of the used software
rule env_export:
    output:
        os.path.join(config["result_path"],'envs','unsupervised_analysis','{env}.yaml'),
    conda:
        "../envs/{env}.yaml"
    resources:
        mem_mb=config.get("mem_small", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_{env}.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        conda env export > {output}
        """