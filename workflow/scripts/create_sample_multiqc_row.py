import os
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def nt_to_aa(sequence):
    return str(Seq(sequence).translate())

def align_seqs(refseq, qseq, both=False):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        temp_input.write(f">ref\n{refseq}\n")
        temp_input.write(f">query\n{qseq}\n")
        temp_input.flush()

        temp_output_name = temp_input.name + ".aln"
        subprocess.run(["clustalw2", "-infile=" + temp_input.name, "-outfile=" + temp_output_name, "-output=fasta", "-type=Protein"], check=True)
    
    for record in SeqIO.parse(temp_output_name, 'fasta'):
        if record.id == 'query':
            queryseq = str(record.seq)
        if record.id == 'ref':
            refseq = str(record.seq)
    if both:
        return queryseq, refseq
    return queryseq

tree_input_path = str(snakemake.input.tree_input)
metadata_path = str(snakemake.input.metadata)
consensus_path = str(snakemake.input.consensus)
complete_path = str(snakemake.input.complete_seq)
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

sequences = {record.id: str(record.seq) for record in SeqIO.parse(complete_path, 'fasta')}

ref_sequences = {}
protein_df = pd.read_excel(metadata_xls, sheet_name='protein_seqs')
for idx, row in protein_df.iterrows():
    gene = row['Gene']
    sequence = row['Sequence']
    org = row['Organism']
    if (gene, org) not in ref_sequences:
        ref_sequences[(gene, org)] = []
    ref_sequences[(gene, org)].append(sequence)

metadata_df = pd.read_excel(metadata_xls, sheet_name='metadata')
known_mutations = [
    {
        'meta_gene': row['Gene'],
        'meta_org': row['Ref. Organism'],
        'meta_pst_pos': row['Ref.Pst aa position'],
        'meta_oldaa': row['WT amino acid'],
        'meta_newaa': row['MUT amino acid'],
        'meta_aapos': int(row['Ref. aa position']),
        'meta_ref': row['Reference']
    }
    for idx, row in metadata_df.iterrows()
]
         
matching_rows = []

for gene, org in ref_sequences:
    for item in known_mutations:
        if gene == item['meta_gene'] and org == item['meta_org']:
            for geneID, sequence in sequences.items():
                if gene == geneID:
                    for aa_refseq in ref_sequences[(gene, org)]:
                        try:
                            aa_seq = nt_to_aa(sequence)
                            aligned_aa_seq = align_seqs(aa_refseq, aa_seq)
                            # If the sequence has all ambiguous aa's, then skip + show didn't amplify note
                            if all(aa == 'X' for aa in aa_seq):
                                matching_rows.append({
                                    'Sample Name': str(snakemake.wildcards.sample),
                                    'Gene': gene,
                                    'Organism': snakemake.wildcards.organism.title(),
                                    'Notes': 'Did not amplify'
                                })
                                raise StopIteration
                        except StopIteration:
                            continue
                        
                        if aligned_aa_seq[item['meta_aapos']-1] in item['meta_newaa']:
                            new_aa = aligned_aa_seq[item['meta_aapos']-1]
                            item['meta_newaa'] = new_aa

                            unaligned_pos = 1
                            aligned_pos = item['meta_aapos'] - 1

                            for i in range(aligned_pos):
                                if aligned_aa_seq[i] != '-':
                                    unaligned_pos += 1

                            new_row = {
                                'Sample Name': str(snakemake.wildcards.sample),
                                'Gene': gene,
                                'Ref.Pst aa position': item['meta_pst_pos'],
                                'WT amino acid': item['meta_oldaa'],
                                'MUT amino acid': new_aa,
                                'Mutation type': 'Non-synonymous',
                                'Ref. Organism': org,
                                'Ref. aa position': item['meta_aapos'],
                                f'{snakemake.wildcards.organism.title()} amino acid position': unaligned_pos,
                                'Reference': item['meta_ref'],
                                'Organism': snakemake.wildcards.organism.title()
                            }
                            if new_row not in matching_rows:
                                matching_rows.append(new_row)

            # Not sure about this, but I'm going to align the {organism}.fna sequence to the samples' translated sequence,
            # and identify any additional mutations against our reference genome.
            # To conserve the mutations already identified above, i'll be aligning against the already aligned sequence of
            # the sample.
            for record in SeqIO.parse(ref_path, 'fasta'):
                if record.id == gene:
                    ref_seq = str(record.seq)
                    aa_ref_seq = nt_to_aa(ref_seq)
                    aligned_ref_aa_seq, aligned_aa_seq = align_seqs(aligned_aa_seq, aa_ref_seq, both=True)
                    for i in range(1, len(aligned_ref_aa_seq) + 1):
                        if aligned_aa_seq[i-1] != aligned_ref_aa_seq[i-1] and aligned_aa_seq[i-1] != 'X':
                            new_aa = aligned_aa_seq[i-1]

                            unaligned_ref_pos = 1
                            unaligned_pos = 1
                            aligned_pos = i - 1
                            for j in range(aligned_pos):
                                if aligned_ref_aa_seq[j] != '-':
                                    unaligned_ref_pos += 1
                                if aligned_aa_seq[j] != '-':
                                    unaligned_pos += 1
                                    
                            new_row = {
                                'Sample Name': str(snakemake.wildcards.sample),
                                'Gene': gene,
                                'WT amino acid': aligned_ref_aa_seq[i-1],
                                'MUT amino acid': new_aa,
                                'Mutation type': 'Non-synonymous',
                                'Ref. Organism': org,
                                'Ref. aa position': unaligned_ref_pos,
                                f'{snakemake.wildcards.organism.title()} amino acid position': unaligned_pos,
                                'Reference': 'N/A',
                                'Organism': snakemake.wildcards.organism.title()
                            }
                            if new_row not in matching_rows:
                                matching_rows.append(new_row) 

matching_df = pd.DataFrame(matching_rows)
matching_df.to_csv(fung_output_path, index=None)
