import gzip
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord

reference_seq = str(snakemake.input.existing_samples)
new_seqs = list(snakemake.input.new_samples)

# Unzip reference_seq and parse it
with gzip.open(reference_seq, 'rt') as ref_file:
    ref_records = {r.id: r for r in SeqIO.parse(ref_file, 'fasta')}

new_records = {}
for new_seq in new_seqs:
    for record in SeqIO.parse(new_seq, 'fasta'):
        tree_bases = 100 * (1 - record.seq.count('N') / len(record.seq))
        if tree_bases > 40:
            new_records[record.id] = record
        else:
            print(f"Skipping {record.id} due to poor amplification")

# new_records = {r.id: r for new_seq in new_seqs for r in SeqIO.parse(new_seq, 'fasta')}

# Concatenate all sequences and align
all_records = {**ref_records, **new_records}
SeqIO.write(all_records.values(), snakemake.output.all_merged, 'fasta')

alignment = AlignIO.read(snakemake.output.all_merged, 'fasta')
common = []

# Identify common positions
for i in range(len(alignment[0])):
    if all((s[i] != '-' and s[i] != 'N') for s in alignment):
        common.append(i)

# Extract common sequences
common_seqs = []
for record in alignment:
    common_seq = ''.join(record.seq[i] for i in common)
    common_seqs.append(SeqRecord(common_seq, id=record.id, description=""))

SeqIO.write(common_seqs, snakemake.output.common_merged, 'fasta')