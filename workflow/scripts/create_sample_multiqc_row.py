import pandas as pd
from Bio import SeqIO
import os
import multiqc
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

def position_range(position):
    if '-' in str(position):
        start, end = map(int, position.split('-'))
        return list(range(start, end + 1))
    else:
        return [int(position)]
    
tree_input_path = str(snakemake.input.tree_input)
metadata_path = str(snakemake.input.metadata)
consensus_path = str(snakemake.input.consensus)
fung_output_path = str(snakemake.output.multiqc_fung)

metadata_df = pd.read_excel(metadata_path)
tree_input_seq = next(str(r.seq) for r in SeqIO.parse(tree_input_path, 'fasta'))
pct_tree_bases = 100 * (1 - tree_input_seq.count('N') / len(tree_input_seq))

groupby_summaries = [pd.DataFrame(
    {'Sample Name': snakemake.wildcards.sample, 'Tree Bases': [pct_tree_bases]}
)]

groupby_cols = ['pool', 'plate']

for i, path in enumerate(snakemake.input.primer_performance_summaries):
    path = str(path)
    groupby_col = os.path.splitext(path.rsplit('_by_')[1])[0]

    primers = pd.read_csv(path)
    if i == 0:
        groupby_summaries.append(
            pd.DataFrame({'All Amplicons': [100 * primers.amplicon_n_used.sum() / primers.amplicon_length.sum()]})
        )
    groupby_summaries.append(pd.DataFrame({
        f'{groupby_col.title()} {val}': [primers.set_index(groupby_col).loc[val]['amplicon_pct_used']]
        for val in primers[groupby_col].values
    }))
pd.concat(groupby_summaries, axis=1).to_csv(snakemake.output.multiqc_row, index=None)

# Adding fungicide resistance genes records to the multiqc report

sequences = {record.id: str(record.seq) for record in SeqIO.parse(consensus_path, 'fasta')}

matching_rows = []
known_mutations = set(metadata_df['Nucleotide Position'])

for idx, row in metadata_df.iterrows():
    gene = row['Gene']
    nuc_pos = position_range(row['Nucleotide Position'])
    old_nuc = row['Old Nucleotide']
    new_nuc = row['New Nucleotide']
    
    for sample_name, sequence in sequences.items():
        if gene == sample_name:
            for pos in nuc_pos:
                print(f'Checking {sample_name} at position {pos} in {gene} (with length {len(sequence)})')
                if sequence[pos - 1] == new_nuc:
                    matching_rows.append({
                        'Sample Name': snakemake.wildcards.sample,
                        **row.to_dict()
                    })

# Check for SNPs not in the metadata
for sample_name, sequence in sequences.items():
    if gene == sample_name:
        for i, nuc in enumerate(sequence):
            if nuc != tree_input_seq[i] and i not in known_mutations:
                matching_rows.append({
                    'Sample Name': snakemake.wildcards.sample,
                    'Gene': sample_name,
                    'Nucleotide Position': i,
                    'Old Nucleotide': tree_input_seq[i],
                    'New Nucleotide': nuc,
                    'Amino Acid Position': '?',
                    'Old AminoAcid': '...',
                    'New AminoAcid': '...'
                })

matching_df = pd.DataFrame(matching_rows)
matching_df.to_csv(fung_output_path, index=None)

