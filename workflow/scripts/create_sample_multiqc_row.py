import os
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def nt_to_aa(sequence):
    return str(Seq(sequence).translate())

def align_seqs(refseq, qseq):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        temp_input.write(f">ref\n{refseq}\n")
        temp_input.write(f">query\n{qseq}\n")
        temp_input.flush()

        temp_output_name = temp_input.name + ".aln"
        subprocess.run(["clustalw2", "-infile=" + temp_input.name, "-outfile=" + temp_output_name, "-output=fasta", "-type=Protein"], check=True)
    
    for record in SeqIO.parse(temp_output_name, 'fasta'):
        if record.id == 'query':
            return str(record.seq)

tree_input_path = str(snakemake.input.tree_input)
metadata_path = str(snakemake.input.metadata)
consensus_path = str(snakemake.input.consensus)
ref_path = str(snakemake.input.reference)
fung_output_path = str(snakemake.output.multiqc_fung)

metadata_xls = pd.ExcelFile(metadata_path)
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

ref_sequences = {}
protein_df = pd.read_excel(metadata_xls, sheet_name='protein_seqs')
for idx, row in protein_df.iterrows():
    gene = row['Gene']
    sequence = row['Sequence']
    if gene in ref_sequences:
        ref_sequences[gene].append(sequence)
    else:
        ref_sequences[gene] = [sequence]

metadata_df = pd.read_excel(metadata_xls, sheet_name='metadata')
matching_rows = []
known_mutations = set(metadata_df['Ref. aa position'])

try:
    for idx, row in metadata_df.iterrows():
        gene = row['Gene']
        ref_org = row['Ref. Organism']
        old_aa = row['WT amino acid']
        new_aa = row['MUT amino acid']
        aa_pos = int(row['Ref. aa position'])
        
        for geneID, sequence in sequences.items():
            if gene == geneID and gene in ref_sequences:
                for aa_refseq in ref_sequences[gene]:
                    aa_seq = nt_to_aa(sequence)
                    aligned_aa_seq = align_seqs(aa_refseq, aa_seq)
                    if all(aa == 'X' for aa in aa_seq):
                        matching_rows.append({
                            'Sample Name': str(snakemake.wildcards.sample),
                            'Gene': gene,
                            'Organism': snakemake.wildcards.organism.title(),
                            'Notes': 'Did not amplify'
                        })
                        raise StopIteration

                    elif aligned_aa_seq[aa_pos-1] in new_aa:
                        new_aa = aligned_aa_seq[aa_pos-1]
                        
                        row['New amino acid'] = new_aa
                        
                        # Get the unaligned position for the P[g/s]t amino acid position
                        unaligned_pos = 0
                        aligned_pos = aa_pos - 1
                        for i in range(1, aligned_pos-1):
                            if aligned_aa_seq[i] != '-':
                                unaligned_pos += 1
                        
                        row[f'{snakemake.wildcards.organism.title()} amino acid position'] = unaligned_pos
                        row[f'Organism'] = snakemake.wildcards.organism.title()
                        new_row = {
                            'Sample Name': str(snakemake.wildcards.sample),
                            **row.to_dict()
                        }
                        if new_row not in matching_rows:
                            matching_rows.append(new_row)
except StopIteration:
    pass

matching_df = pd.DataFrame(matching_rows)
matching_df.to_csv(fung_output_path, index=None)
