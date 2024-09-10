import os
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# snakemake_organism = snakemake.wildcards.organism
snakemake_organism = 'pgt'

def nt_to_aa(sequence):
    return str(Seq(sequence).translate())

def align_seqs(refseq, qseq, both=True):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        temp_input.write(f">ref\n{refseq}\n>query\n{qseq}\n")
        temp_input.flush()
        temp_output_name = temp_input.name + ".aln"
        # Need to change to work on Mac -- can use MAFFT instead but need to handle the output differently 
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

# Adding fungicide resistance genes records to the multiqc report

metadata_xls = pd.ExcelFile(str(snakemake.input.metadata))

sequences = {record.id: str(record.seq) for record in SeqIO.parse(str(snakemake.input.consensus), 'fasta')}
ref_sequences = {(row['Gene'], row['Organism']): row['Sequence'] for idx, row in load_data(metadata_xls, 'protein_seqs').iterrows()}
known_mutations = [row.to_dict() for idx, row in load_data(metadata_xls, 'metadata').iterrows()]
variants = [row.to_dict() for idx, row in load_data(metadata_xls, 'variants').iterrows()]

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
                        new_row = {'Sample Name': str(snakemake.wildcards.sample), 'Gene': gene, 'Organism': snakemake_organism.title(), 'Notes': 'Did not amplify'}
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
                            'WT amino acid': item['WT amino acid'], 'MUT amino acid': new_aa,
                            'Ref. Organism': org, 'Ref. aa position': item['Ref. aa position'], 
                            f'Pgt amino acid position': unaligned_sample_pos, 'Reference': item['Reference'],
                            'Organism': snakemake_organism.title()
                        }
                        if new_row not in matching_rows:
                            matching_rows.append(new_row)
                        
    ref_seq = next(str(record.seq) for record in SeqIO.parse(str(snakemake.input.reference), 'fasta') if record.id == gene)
    aa_ref_seq = nt_to_aa(ref_seq)
    aa_seq = nt_to_aa(sequences[gene])
    aligned_aa_seq, aligned_ref_aa_seq = align_seqs(aa_ref_seq, aa_seq)

    ref_alignment_map = alignment_mapping(aligned_ref_aa_seq)
    sample_alignment_map = alignment_mapping(aligned_aa_seq)

    for i, (ref_aa, seq_aa) in enumerate(zip(aligned_ref_aa_seq, aligned_aa_seq), start=1):
        if (ref_aa != seq_aa) and seq_aa != 'X': # we're also doing insertions/deletions
            unaligned_ref_pos = ref_alignment_map[i - 1]
            unaligned_pos = sample_alignment_map[i - 1]
            
            # Check if the mutation is a different variant of the same gene
            variant_found = False
            for item in variants:
                if item['Gene'] == gene and unaligned_ref_pos == item['Posn'] and \
                    ((seq_aa == item['Var1']) or (seq_aa == item['Var2']) or (seq_aa == item['Var3'])):
                    variant_found = True
            if variant_found:
                continue

            new_row = {
                'Sample Name': str(snakemake.wildcards.sample), 'Gene': gene, 'WT amino acid': ref_aa, 'MUT amino acid': seq_aa, 
                'Ref. Organism': f'{snakemake_organism.title()}', 'Ref. aa position': unaligned_ref_pos,
                f'{snakemake_organism.title()} amino acid position': unaligned_pos, 'Reference': 'N/A', 
                'Organism': snakemake_organism.title()
            }
            if new_row not in matching_rows:
                matching_rows.append(new_row)

pd.DataFrame(matching_rows).to_csv(snakemake.output.multiqc_fung, index=None)

# Adding R-genes and Avr records to the multiqc report

def out_vir(sequence):
    amb_count = sequence.count('N')
    pct_aligned = (len(sequence) - amb_count)*100 / len(sequence)
    print(pct_aligned)
    if pct_aligned > 50:
        return {'Sample Name': snakemake.wildcards.sample,'Gene ID': geneID, 'Percentage Aligned': pct_aligned}
    return None

avr_rows = []
sr_rows = []

for geneID, sequence in sequences.items():
    if geneID.startswith('Avr'):
        print(geneID)
        new_row = out_vir(sequence)
        if new_row:
            avr_rows.append(new_row)
    elif geneID.startswith('Sr'):
        print(geneID)
        new_row = out_vir(sequence)
        if new_row:
            sr_rows.append(new_row)

pd.DataFrame(avr_rows).to_csv(snakemake.output.multiqc_avr, index=None)
pd.DataFrame(sr_rows).to_csv(snakemake.output.multiqc_rgene, index=None)
