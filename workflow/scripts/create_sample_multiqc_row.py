import os
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def nt_to_aa(sequence):
    return str(Seq(sequence).translate())

def align_seqs(refseq, qseq, both=True):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        temp_input.write(f">ref\n{refseq}\n>query\n{qseq}\n")
        temp_input.flush()
        temp_output_name = temp_input.name + ".aln"
        subprocess.run(["clustalw2", "-infile=" + temp_input.name, "-outfile=" + temp_output_name, "-output=fasta", "-type=Protein", '-quiet=stdout'], check=True)
    
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(temp_output_name, 'fasta')}
    return (sequences['query'], sequences['ref']) if both else sequences['query']

def load_data(file_path, sheet_name=None):
    return pd.read_excel(file_path, sheet_name=sheet_name) if sheet_name else pd.read_csv(file_path)

def alignment_mapping(aligned_seq):
    alignment_map = {}
    unaligned_pos = 0
    for aligned_pos, aa in enumerate(aligned_seq):
        if aa != '-':
            unaligned_pos += 1
        alignment_map[aligned_pos] = unaligned_pos
    return alignment_map

tree_input_seq = next(str(r.seq) for r in SeqIO.parse(str(snakemake.input.tree_input), 'fasta'))
pct_tree_bases = 100 * (1 - tree_input_seq.count('N') / len(tree_input_seq))

metadata_xls = pd.ExcelFile(str(snakemake.input.metadata))
groupby_summaries = [pd.DataFrame({'Sample Name': snakemake.wildcards.sample, 'Tree Bases': [pct_tree_bases]})]

for i, path in enumerate(snakemake.input.primer_performance_summaries):
    primers = load_data(str(path))
    groupby_col = os.path.splitext(str(path).rsplit('_by_')[1])[0]
    if i == 0:
        groupby_summaries.append(pd.DataFrame({'All Amplicons': [100 * primers.amplicon_n_used.sum() / primers.amplicon_length.sum()]}))
    groupby_summaries.append(pd.DataFrame({f'{groupby_col.title()} {val}': [primers.set_index(groupby_col).loc[val]['amplicon_pct_used']] for val in primers[groupby_col].values}))

pd.concat(groupby_summaries, axis=1).to_csv(snakemake.output.multiqc_row, index=None)

# Adding fungicide resistance genes records to the multiqc report

sequences = {record.id: str(record.seq) for record in SeqIO.parse(str(snakemake.input.complete_seq), 'fasta')}
ref_sequences = {(row['Gene'], row['Organism']): row['Sequence'] for idx, row in load_data(metadata_xls, 'protein_seqs').iterrows()}
known_mutations = [row.to_dict() for idx, row in load_data(metadata_xls, 'metadata').iterrows()]

matching_rows = []

for gene, org in ref_sequences:
    aa_refseq = ref_sequences[(gene, org)]
    for item in known_mutations:
        if gene == item['Gene'] and org == item['Ref. Organism']:
            for geneID, sequence in sequences.items():
                if gene == geneID:
                    aa_seq = nt_to_aa(sequence)
                    aligned_aa_seq, aligned_ref_aa_seq = align_seqs(aa_refseq, aa_seq)

                    ref_alignment_map = alignment_mapping(aligned_ref_aa_seq)
                    sample_alignment_map = alignment_mapping(aligned_aa_seq)
                    
                    if all(aa == 'X' for aa in aa_seq):
                        new_row = {'Sample Name': str(snakemake.wildcards.sample), 'Gene': gene, 'Organism': snakemake.wildcards.organism.title(), 'Notes': 'Did not amplify'}
                        if new_row not in matching_rows:
                            matching_rows.append(new_row)
                        break

                    ref_pos = item['Ref. aa position'] - 1
                    aligned_ref_pos = next((pos for pos, unaligned_pos in ref_alignment_map.items() if unaligned_pos == ref_pos + 1), None)
                    
                    if aligned_ref_pos is not None and aligned_aa_seq[aligned_ref_pos] in item['MUT amino acid'] and aligned_ref_aa_seq[aligned_ref_pos] != '-':
                        new_aa = aligned_aa_seq[aligned_ref_pos]
                        unaligned_sample_pos = sample_alignment_map[aligned_ref_pos]
                        new_row = {
                            'Sample Name': str(snakemake.wildcards.sample), 'Gene': gene, 'Ref.Pst aa position': item['Ref.Pst aa position'], 
                            'WT amino acid': item['WT amino acid'], 'MUT amino acid': new_aa, 'Mutation type': 'Non-synonymous',
                            'Ref. Organism': org, 'Ref. aa position': item['Ref. aa position'], 
                            f'{snakemake.wildcards.organism.title()} amino acid position': unaligned_sample_pos, 'Reference': item['Reference'],
                            'Organism': snakemake.wildcards.organism.title()
                        }
                        if new_row not in matching_rows:
                            matching_rows.append(new_row)

    # Can also check against the consensus sequence of the P[g/s]t reference, but initial tests showed great differences between the complete sequence and the references, so might be better to leave out
    
    # ref_seq = next(str(record.seq) for record in SeqIO.parse(str(snakemake.input.reference), 'fasta') if record.id == gene)
    # aa_ref_seq = nt_to_aa(ref_seq)
    # aa_seq = nt_to_aa(sequences[gene])
    # aligned_aa_seq, aligned_ref_aa_seq = align_seqs(aa_ref_seq, aa_seq)

    # ref_alignment_map = alignment_mapping(aligned_ref_aa_seq)
    # sample_alignment_map = alignment_mapping(aligned_aa_seq)

    # for i, (ref_aa, seq_aa) in enumerate(zip(aligned_ref_aa_seq, aligned_aa_seq), start=1):
    #     if (ref_aa != seq_aa) and seq_aa != 'X': # we're also doing insertions/deletions
    #         unaligned_ref_pos = ref_alignment_map[i - 1]
    #         unaligned_pos = sample_alignment_map[i - 1]
    #         new_row = {
    #             'Sample Name': str(snakemake.wildcards.sample), 'Gene': gene, 'WT amino acid': ref_aa, 'MUT amino acid': seq_aa, 
    #             'Mutation type': 'Non-synonymous', 'Ref. Organism': f'{snakemake.wildcards.organism.title()}', 'Ref. aa position': unaligned_ref_pos,
    #             f'{snakemake.wildcards.organism.title()} amino acid position': unaligned_pos, 'Reference': 'N/A', 
    #             'Organism': snakemake.wildcards.organism.title()
    #         }
    #         if new_row not in matching_rows:
    #             matching_rows.append(new_row)

pd.DataFrame(matching_rows).to_csv(snakemake.output.multiqc_fung, index=None)