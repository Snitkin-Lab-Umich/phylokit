# Authors: Kyle Gontjes & Dhatri Badri
# Date: 12-20-2024

configfile: "config/config.yaml"

import pandas as pd
import os
import re

PREFIX = config["prefix"]

if not os.path.exists("results/" + PREFIX):
    try:
        os.makedirs("results/" + PREFIX)
    except OSError as e:
        print(f"Error creating directory: {e}")

# Rule all information
rule all:
    input:

rule gubbins:
    input: 

    output:
        recombination_predictions_gff = "results/{{prefix}}/gubbins/{{prefix}}.recombination_predictions.gff",
        recombination_predictions_embl = "results/{{prefix}}/gubbins/{{prefix}}.recombination_predictions.embl",
        base_reconstruction = "results/{{prefix}}/gubbins/{{prefix}}.branch_base_reconstruction.embl",
        snp_mutation_summary = "results/{{prefix}}/gubbins/{{prefix}}.summary_of_snp_distribution.vcf",
        per_branch_stats = "results/{{prefix}}/gubbins/{{prefix}}.per_branch_statistics.csv",
        snp_mutation_summary = "results/{{prefix}}/gubbins/{{prefix}}.filtered_polymorphic_sites.fasta",
        filtered_polymorphic_sites_phylip_fmt = "results/{{prefix}}/gubbins/{{prefix}}.filtered_polymorphic_sites.phylip",
        newick_tree = "results/{{prefix}}/gubbins/{{prefix}}.final_tree.tre",
        phylo_tree_node_labels = "results/{{prefix}}/gubbins/{{prefix}}.node_labelled.final_tree.tre",
        log = "results/{{prefix}}/gubbins/{{prefix}}.log"
    params:
        prefix = config["prefix"],
        threads = config["threads"],
        first_tree_algorithm = config["first_tree_algorithm"],
        tree_builder = config["tree_builder"],
        tree_builder_model = config["tree_builder_model"],
        outgroup = config["outgroup"],
        alignment = config["alignment"],
        
    singularity:
        "docker://staphb/gubbins:3.3.5"
    shell:
        "run_gubbins.py --prefix {params.prefix} --threads {params.threads} --verbose --no-cleanup --first-tree-builder {params.first_tree_algorithm} --first-model {params.tree_builder_model} --tree-builder {params.tree_builder_model} --outgroup {params.outgroup} {params.alignment}"

rule mask_gubbins:
    input:
        recombination_predictions_gff = lambda wildcards: expand(f"results/{wildcards.prefix}/gubbins/{wildcards.prefix}.recombination_predictions.gff"),
    output:
        gubbins_masked = "results/{{prefix}}/phylo_analysis/gubbins/{{ref_name}}_gubbins_masked.fa"
        masked_recomb_pos = "results/{{prefix}}/phylo_analysis/genome_alignment/{{ref_name}}_masked_recomb_positions.txt"
    params:
        prefix =  config["prefix"],
        outdir = "results/{prefix}/gubbins"
    wrapper:
        "file:wrapper_functions/mask_gubbins"

rule iqtree:
    input:
        gubbins_masked = lambda wildcards: expand(f"results/{wildcards.prefix}/gubbins/{wildcards.prefix}.filtered_polymorphic_sites.fasta"),
    output:
        ml_tree_newick_fmt = "results/{{prefix}}/IQtree/{{prefix}}.treefile"
        iqtree_report = "results/{{prefix}}/IQtree/.iqtree"
    params:
       iqtree_model = config["iqtree_model"],
       bootstrap_count = config["bootstrap_count"],
       outgroup = config["outgroup"],
       num_unsuccessful_iterations = config["num_unsuccessful_iterations"],
       prefix = config["prefix"]
    singularity:
        "docker://staphb/iqtree2:2.3.6"
    shell:
        "iqtree -s {input.gubbins_masked} -nt AUTO -m {params.iqtree_model} -B {params.bootstrap_count} -bnni -o {params.outgroup} -nstop {params.num_unsuccessful_iterations} -pre {params.prefix}" 