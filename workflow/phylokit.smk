# Authors: Kyle Gontjes & Dhatri Badri
# Date: 04-30-2025

configfile: "config/config.yaml"

import pandas as pd
import os
import re
 
PREFIX = config["prefix"]
RUN = config["run_name"]

if not os.path.exists("results/"):
    try:
        os.makedirs("results/")
    except OSError as e:
        print(f"Error creating directory: {e}")

# Rule all information
rule all:
    input:
        ml_tree_newick_fmt = expand("results/{run}/IQtree/{prefix}.treefile", run=RUN, prefix=PREFIX)

rule gubbins:
    input: 
        alignment = config["alignment"]
    output:
        recombination_predictions_gff = "results/{run}/gubbins/{prefix}.recombination_predictions.gff"
    params:
        pwd = "results/{run}/gubbins/", 
        file_prefix = config["prefix"],
        threads = config["threads"],
        first_tree_algorithm = config["first_tree_algorithm"],
        first_model = config["first_model"], 
        outgroup = config["outgroup"] 
    singularity:
        "docker://staphb/gubbins:3.3.5"
    shell:
        """
        cd {params.pwd} && \
        run_gubbins.py --prefix {params.file_prefix} --threads {params.threads} --verbose --first-tree-builder {params.first_tree_algorithm} --first-model {params.first_model} --best-model --outgroup {params.outgroup} {input.alignment}
        """

rule mask_gubbins:
    input:
        recombination_predictions_gff = lambda wildcards: expand(f"results/{wildcards.run}/gubbins/{wildcards.prefix}.recombination_predictions.gff")
    output:
        gubbins_masked_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_gubbins_masked_var_sites.fa"
    params:
        outdir = "results/{run}/gubbins_masked",
        prefix =  config["prefix"],
        alignment = config["alignment"]
    wrapper:
        "file:wrapper_functions/mask_gubbins" 

rule iqtree:
    input:
        gubbins_masked_var_sites = lambda wildcards: expand(f"results/{wildcards.run}/gubbins_masked/{wildcards.prefix}_gubbins_masked_var_sites.fa") 
    output:
        ml_tree_newick_fmt = f"results/{{run}}/IQtree/{{prefix}}.treefile" 
    params:
       iqtree_model = config["iqtree_model"],
       bootstrap_count = config["bootstrap_count"],
       outgroup = config["outgroup"],
       num_unsuccessful_iterations = config["num_unsuccessful_iterations"],
       iqtree_prefix = expand(str(config["prefix"])),
       iqtree_dir = "results/{run}/IQtree"
    conda:
        "envs/iqtree.yaml"
    log:
        iqtree_log = "logs/{run}/IQtree/{prefix}.iqtree.log"
    shell:
        "iqtree -s {input.gubbins_masked_var_sites} -nt AUTO -m {params.iqtree_model} -B {params.bootstrap_count} -o {params.outgroup} -nstop {params.num_unsuccessful_iterations} -pre {params.iqtree_dir}/{params.iqtree_prefix} &> {log.iqtree_log}" 
