# conda environment export rule to document the exact versions and builds of the used software
rule env_export:
    output:
        report(os.path.join(result_path,'envs','{env}.yaml'),
                      caption="../report/software.rst", 
                      category="Software", 
                      subcategory="{}_{}".format(config["project_name"], module_name)
                     ),
    conda:
        "../envs/{env}.yaml"
    resources:
        mem_mb=1000,
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_{env}.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        conda env export > {output}
        """
        
# export used config file for documentation and reproducibility
rule config_export:
    output:
        configs = report(os.path.join(result_path,'configs','{}_config.yaml'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_{}".format(config["project_name"], module_name)
                        )
    resources:
        mem_mb=1000,
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","config_export.log"),
    params:
        partition=config.get("partition"),
    run:
        with open(output["configs"], 'w') as outfile:
            yaml.dump(config, outfile)

# export used annotation file for documentation and reproducibility         
rule annot_export:
    input:
        config["annotation"],
    output:
        annot = report(os.path.join(result_path,'configs','{}_annot.csv'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_{}".format(config["project_name"], module_name)
                        )
    resources:
        mem_mb=1000,
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","annot_export.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        cp {input} {output}
        """
