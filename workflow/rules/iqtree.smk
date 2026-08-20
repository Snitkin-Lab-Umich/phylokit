rule iqtree:
    input:
        gubbins_masked_var_sites = lambda wildcards: expand(f"results/{wildcards.run}/gubbins_masked/{wildcards.prefix}_gubbins_masked_var_sites.fa") 
    output:
        ml_tree_newick_fmt = f"results/{{run}}/iqtree/{{prefix}}.treefile" 
    params:
        iqtree_model=config["iqtree_model"],
        bootstrap_count=config["bootstrap_count"],
        num_unsuccessful_iterations=config["num_unsuccessful_iterations"],
        iqtree_prefix=lambda wildcards, output: os.path.splitext(output.ml_tree_newick_fmt)[0],
        outgroup_flag=lambda wildcards: (
            f"-o {config['outgroup']}" if config.get("outgroup") else ""
        ),
        bnni_flag=lambda wildcards: "-bnni" if config.get("bnni", True) else "",
        mem_gb=config["mem_mb_iqtree"] // 1000,
    conda:
        "../envs/iqtree.yaml"
    log:
        "logs/{run}/iqtree/{prefix}.iqtree.log"
    benchmark: 
        "benchmarks/{run}/iqtree/{prefix}.benchmark.tsv"
    resources:
        mem_mb = config["mem_mb_iqtree"],
        runtime = config["runtime_iqtree"]
    threads:
        config["threads"]
    shell:
        "iqtree -s {input.gubbins_masked_var_sites} -st DNA -T {threads} -m {params.iqtree_model} -mem {params.mem_gb}G -bb {params.bootstrap_count} -nstop {params.num_unsuccessful_iterations} -pre {params.iqtree_prefix} {params.outgroup_flag} {params.bnni_flag} &> {log}"