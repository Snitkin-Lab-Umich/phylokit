# phylokit
Snakemake pipeline to generate a phylogenetic tree from an alignment 

# Installation
```
git clone https://github.com/Snitkin-Lab-Umich/phylokit.git
```

# Simple command
```
module load snakemake

module load singularity

module load mamba

snakemake -s phylokit.smk --use-conda  --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}  --output=slurm_out/slurm-%j.out" --conda-frontend mamba --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 60
```