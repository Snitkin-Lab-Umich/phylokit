from snakemake.shell import shell
import os 
import sys
import pandas as pd 
from collections import defaultdict
import Bio
from Bio import AlignIO, SeqIO  
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

outdir = snakemake.params[0]
prefix = snakemake.params[1]
alignment = snakemake.params[2] 
gff_file = snakemake.input[0]

def is_gff_effectively_empty(gff_file):
    if os.stat(gff_file).st_size == 0:
        return True
    with open(gff_file) as f:
        for line in f:
            if not line.startswith("#"):
                return False
    return True

def mask_recombinant_regions(gff_file, alignment, outdir, prefix):
    # read in alignment and gubbins gff file
    gff = pd.read_csv(gff_file, sep='\t', skiprows=2, header=None)
    aln = AlignIO.read(alignment, 'fasta')

    # get indices/positions of recombinant regions
    recomb_regions = defaultdict(list)
    for row in gff.iterrows():
        start = row[1][3]
        end = row[1][4]
        region = list(range(start, end + 1))
        taxa = row[1][8].split(';')[2]
        taxa = taxa.replace('taxa="', '').replace('"', '')
        taxa = list(taxa.split())
        for isolate in taxa:
            for position in region:
                recomb_regions[isolate].append(position)

    # mask recombinant regions
    sample_masked_indices = defaultdict(list)
    new_aln = []
    for record in aln:
        seq_str = list(str(record.seq))
        masked_indices = recomb_regions.get(record.id, [])
        for index in masked_indices:
            seq_str[index - 1] = 'N'
        seq_str = ''.join(seq_str)
        new_record = SeqRecord(Seq(seq_str), id=record.id, description="")
        sample_masked_indices[record.id] = masked_indices
        new_aln.append(new_record)

    # write new FASTA file with recombinant regions masked
    fasta_outfile = os.path.join(outdir, prefix + "_gubbins_masked.fa")
    text_outfile = os.path.join(outdir, prefix + "_masked_recomb_positions.txt")
    var_site_outfile = os.path.join(outdir, prefix + "_gubbins_masked_var_sites.fa")

    with open(fasta_outfile, 'w') as handle:
        SeqIO.write(new_aln, handle, 'fasta')

    with open(text_outfile, 'w') as handle:
        for sample, positions in sample_masked_indices.items():
            line = str(sample) + '\t' + ','.join(map(str, positions)) + '\n'
            handle.write(line)

    cmd = f'snp-sites {fasta_outfile} -m -o {var_site_outfile}'
    os.system(cmd)

def main(gff_file, alignment, outdir, prefix):
    var_site_outfile = os.path.join(outdir, prefix + "_gubbins_masked_var_sites.fa")
    if is_gff_effectively_empty(gff_file):
        print("GFF is empty or has no data. Running snp-sites on original alignment.")
        cmd = f'snp-sites {alignment} -m -o {var_site_outfile}'
        os.system(cmd)
    else:
        print("GFF contains data. Masking recombinant regions before phylogenetic anaysis.")
        mask_recombinant_regions(gff_file, alignment, outdir, prefix)

main(gff_file, alignment, outdir, prefix)

