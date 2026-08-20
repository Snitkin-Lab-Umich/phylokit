# Author: Kyle Gontjes + Ali Pirani

rule mask_gubbins:
    input:
        recombination_predictions_gff = "results/{run}/gubbins/{prefix}.recombination_predictions.gff",
    output:
        gubbins_masked_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_gubbins_masked_var_sites.fa",
        pre_gubbins_masked_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_pre_gubbins_masked_var_sites.fa",
        gubbins_masked_core_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_gubbins_masked_core_var_sites.fa",
        pre_gubbins_masked_core_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_pre_gubbins_masked_core_var_sites.fa",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.gubbins_masked_var_sites), 
        prefix =  config["prefix"],
        alignment = config["alignment"]
    resources:
        mem_mb=config["mem_mb_gubbins_mask"],
        runtime=config["runtime_gubbins_mask"],
    log:
        "logs/{run}/gubbins_masked/{prefix}.mask_gubbins.log",
    benchmark:
        "benchmarks/{run}/gubbins_masked/{prefix}.benchmark.tsv",
    conda:
        "../envs/mask_gubbins.yaml"
    script:
        "../scripts/mask_gubbins.py"