# Author: Kyle Gontjes + Ali Pirani

rule gubbins:
    input: 
        alignment = config["alignment"]
    output:
        recombination_predictions_gff = "results/{run}/gubbins/{prefix}.recombination_predictions.gff",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.recombination_predictions_gff),
        prefix=lambda wildcards: wildcards.prefix,
        first_tree_algorithm = config["first_tree_algorithm"],
        first_model = config["first_model"], 
        tree_builder = config['tree_builder'],
        filter_percentage = config['filter_percentage'],
        outgroup_flag=lambda wildcards: (
            f"--outgroup {config['outgroup']}" if config.get("outgroup") else ""
        ),
    log:
        "logs/{run}/gubbins/{prefix}.gubbins.log",
    benchmark: 
        "benchmarks/{run}/gubbins/{prefix}.benchmark.tsv"
    resources:
        mem_mb=config["mem_mb_gubbins"],
        runtime=config["runtime_gubbins"],
    threads:
        config["threads"]
    container:
        "docker://staphb/gubbins:3.4.1"
    shell:
        """
        LOGFILE=$(realpath {log})
        mkdir -p {params.outdir}
        cd {params.outdir}
        run_gubbins.py --version > $LOGFILE 2>&1
        run_gubbins.py --prefix {params.prefix} --threads {threads} --verbose --filter-percentage {params.filter_percentage} --first-tree-builder {params.first_tree_algorithm} --first-model {params.first_model} --tree-builder {params.tree_builder} --best-model {params.outgroup_flag} {input.alignment}  >> $LOGFILE 2>&1
        """