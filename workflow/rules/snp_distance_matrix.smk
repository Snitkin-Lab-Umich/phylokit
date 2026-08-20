rule snp_distance_matrix:
    input:
        pre_gubbins_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_pre_gubbins_masked_var_sites.fa",
        pre_gubbins_core_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_pre_gubbins_masked_core_var_sites.fa",
        gubbins_masked_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_gubbins_masked_var_sites.fa",
        gubbins_masked_core_var_sites = f"results/{{run}}/gubbins_masked/{{prefix}}_gubbins_masked_core_var_sites.fa",
    output:
        pre_gubbins_noncore_snp_matrix = "results/{run}/snp_distance_matrices/{prefix}_pre_gubbins_noncore_snp_matrix.csv", 
        pre_gubbins_core_snp_matrix = "results/{run}/snp_distance_matrices/{prefix}_pre_gubbins_core_snp_matrix.csv",
        gubbins_masked_noncore_snp_matrix = "results/{run}/snp_distance_matrices/{prefix}_gubbins_masked_noncore_snp_matrix.csv",
        gubbins_masked_core_snp_matrix = "results/{run}/snp_distance_matrices/{prefix}_gubbins_masked_core_snp_matrix.csv",
    conda:
        "../envs/snp_distance_matrix.yaml"
    log:
        "logs/{run}/snp_distance_matrices/{prefix}.snp_distance_matrix.log"
    benchmark: 
        "benchmarks/{run}/snp_distance_matrices/{prefix}.benchmark.tsv"
    resources:
        mem_mb = config["mem_mb_snp_distance_matrix"],
        runtime = config["runtime_snp_distance_matrix"]
    threads:
        config["threads"] 
    shell:
        """       
        echo "Starting snp_distance_matrix rule $(date)" > {log}
        
        # Pairwise SNP Matrix
        ## Pre-Gubbins filtered
        pairsnp -c -t {threads} {input.pre_gubbins_var_sites} > {output.pre_gubbins_noncore_snp_matrix} 2>> {log}
        ## Gubbins filtered
        pairsnp -c -t {threads} {input.gubbins_masked_var_sites} > {output.gubbins_masked_noncore_snp_matrix} 2>> {log}
        # Core SNP Matrix
        ## Pre-Gubbins filtered
        pairsnp -c -t {threads} {input.pre_gubbins_core_var_sites} > {output.pre_gubbins_core_snp_matrix} 2>> {log}    
        ## Gubbins filtered
        pairsnp -c -t {threads} {input.gubbins_masked_core_var_sites} > {output.gubbins_masked_core_snp_matrix} 2>> {log}       
        
        echo "Ending snp_distance_matrix rule $(date)" >> {log}
        """