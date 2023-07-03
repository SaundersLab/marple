import pandas as pd
from Bio.SeqIO import parse
import sys

ref_path = str(snakemake.input.ref)
snps_paths = snakemake.input.snps
output = str(snakemake.output.snp_counts)

n_snps = pd.DataFrame(index=[r.id for r in parse(ref_path, 'fasta')])
n_snps.index.name = 'seqid'
for snps_path in snps_paths:
    snps_path = str(snps_path)
    df = pd.read_csv(snps_path, sep='\t')
    snps = df[df.ref != df.alleles]
    sample = snps_path.split('/')[-1].split('.')[0]
    sample_n_snps = pd.DataFrame({sample: snps.groupby('seqid').size()})
    n_snps = n_snps.join(sample_n_snps).fillna(0)
for col in n_snps:
    n_snps[col] = n_snps[col].astype(int)
n_snps.reset_index().sort_values(by='seqid').to_csv(output, index=None)
