import os
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

snakemake_organism = snakemake.wildcards.organism

def nt_to_aa(sequence,mtDNA=False):
    if mtDNA:
        return str(Seq(sequence).translate(table=4))
    else:
        return str(Seq(sequence).translate())

def align_seqs(refseq, qseq, both=True):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        temp_input.write(f">ref\n{refseq}\n>query\n{qseq}\n")
        temp_input.flush()
        temp_output_name = temp_input.name + ".aln"
        # Need to change to work on Mac -- can use MAFFT instead but need to handle the output differently 
        subprocess.run(["clustalw2", "-infile=" + temp_input.name, "-outfile=" + temp_output_name, "-output=fasta", "-type=Protein", '-quiet=stdout'], check=True)
    
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(temp_output_name, 'fasta')}
    os.remove(temp_input.name)
    os.remove(temp_output_name)
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

def load_fungicide_gene_sequences(fasta_file, reps=False):
    """Load fungicide gene sequences from FASTA file into a dictionary."""
    sequences = {}

    if reps:
        sample_xls = pd.ExcelFile(str(snakemake.input.sample_metadata))
        representative_map = {row['tree_name']: row['tree_new_name'] for _, row in load_data(sample_xls, "metadata").iterrows()}
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene = record.description.split("(")[1].split("-")[0].replace(" ", "")
            sample = record.description.split("(")[1].split("-")[1].replace(")", "").replace(" ", "")
            
            new_name = representative_map.get(sample, sample)
            
            sequences.setdefault(gene, {})[new_name] = str(record.seq).upper()
        return sequences

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq).upper()

    return sequences

def add_row(rows, row, seen):
    key = tuple(sorted(row.items()))
    if key not in seen:
        rows.append(row)
        seen.add(key)

def get_ref_seq(gene):
    return next(str(record.seq) for record in SeqIO.parse(str(snakemake.input.reference), "fasta") if record.id == gene)

def find_aligned_pos(alignment_map, unaligned_pos):
    return next((pos for pos, u_pos in alignment_map.items() if u_pos == unaligned_pos), None)

def find_representatives_with_mutation(rep_sequences, gene, position, mutation_aa, ref_organism=None):
    """Find which representatives have the mutation at the given position."""
    representatives_with_mut = []
    
    if gene not in rep_sequences:
        return representatives_with_mut
    
    for sample_name, seq_data in rep_sequences[gene].items():
        # Use organism-specific alignment if provided, otherwise use default
        if ref_organism and ref_organism in seq_data.get("alignments_by_organism", {}):
            alignment_data = seq_data["alignments_by_organism"][ref_organism]
        else:
            alignment_data = seq_data
        
        alignment_map = alignment_data["alignment_map"]
        aligned_seq = alignment_data["aligned_seq"]
        
        # Find aligned position from unaligned position
        aligned_pos = find_aligned_pos(alignment_map, position)
        if aligned_pos is not None and aligned_pos < len(aligned_seq):
            if aligned_seq[aligned_pos] == mutation_aa:
                representatives_with_mut.append(sample_name)
    
    return representatives_with_mut

# Adding fungicide resistance genes records to the multiqc report

metadata_xls = pd.ExcelFile(str(snakemake.input.metadata))

sequences = {record.id: str(record.seq) for record in SeqIO.parse(str(snakemake.input.consensus), 'fasta')}
ref_sequences = {(row['Gene'], row['Organism']): row['Sequence'] for idx, row in load_data(metadata_xls, 'protein_seqs').iterrows()}
known_mutations = [row.to_dict() for idx, row in load_data(metadata_xls, 'metadata').iterrows()]
variants = [row.to_dict() for idx, row in load_data(metadata_xls, 'variants').iterrows()]

matching_rows = []
seen_rows = set()

# Load representative sequences
rep_sequences = load_fungicide_gene_sequences(snakemake.input.representatives, reps=True)

# Align sequences against reference, store aligned sequences for each sample and gene
# Also store alignments for each reference organism
for gene in rep_sequences:
    ref_seq = get_ref_seq(gene)
    aa_ref_seq = nt_to_aa(ref_seq)
    for sample, seq in rep_sequences[gene].items():
        aligned_seq, aligned_ref = align_seqs(aa_ref_seq, seq)
        rep_sequences[gene][sample] = {
            "aligned_seq": aligned_seq,
            "aligned_ref": aligned_ref,
            "alignment_map": alignment_mapping(aligned_ref),
            "alignments_by_organism": {}
        }
    
    # Also align representatives against each reference organism's sequence
    for (ref_gene, ref_org), ref_aa_seq in ref_sequences.items():
        if ref_gene == gene:
            for sample, seq_data in rep_sequences[gene].items():
                aligned_seq, aligned_ref = align_seqs(ref_aa_seq, seq_data["aligned_seq"].replace("-", ""))
                seq_data["alignments_by_organism"][ref_org] = {
                    "aligned_seq": aligned_seq,
                    "aligned_ref": aligned_ref,
                    "alignment_map": alignment_mapping(aligned_ref)
                }

for gene, org in ref_sequences:
    aa_refseq = ref_sequences[(gene, org)]
    if gene not in sequences:
        continue
    
    if gene == "Cob":
        aa_seq = nt_to_aa(sequences[gene],mtDNA=True)
    else:
        aa_seq = nt_to_aa(sequences[gene])
    if all(aa == "X" for aa in aa_seq):
        add_row(
            matching_rows,
            {"Sample Name": str(snakemake.wildcards.sample), "Gene": gene, "Organism": snakemake_organism.title(), "Notes": "Did not amplify"},
            seen_rows,
        )
        continue

    aligned_aa_seq, aligned_ref_aa_seq = align_seqs(aa_refseq, aa_seq)
    ref_alignment_map = alignment_mapping(aligned_ref_aa_seq)
    sample_alignment_map = alignment_mapping(aligned_aa_seq)

    for item in known_mutations:
        if gene != item["Gene"] or org != item["Ref. Organism"]:
            continue

        aligned_ref_pos = find_aligned_pos(ref_alignment_map, item["Ref. aa position"])
        if aligned_ref_pos is None:
            continue

        if aligned_ref_aa_seq[aligned_ref_pos] == "-":
            continue

        if aligned_aa_seq[aligned_ref_pos] in item["MUT amino acid"]:
            new_aa = aligned_aa_seq[aligned_ref_pos]
            unaligned_sample_pos = sample_alignment_map[aligned_ref_pos]
            
            # Find which representatives have this mutation against the reference organism
            reps_with_mut = find_representatives_with_mutation(rep_sequences, gene, item["Ref. aa position"], new_aa, ref_organism=org)
            total_reps = len(rep_sequences.get(gene, {}))
            
            add_row(
                matching_rows,
                {
                    "Sample Name": str(snakemake.wildcards.sample),
                    "Gene": gene,
                    "WT amino acid": item["WT amino acid"],
                    "MUT amino acid": new_aa,
                    "Ref. Organism": org,
                    "Ref. aa position": item["Ref. aa position"],
                    "Sample amino acid position": unaligned_sample_pos,
                    "Reference": item["Reference"],
                    "Organism": str(snakemake.wildcards.organism).title(),
                    "Representatives with mutation": "All representatives" if total_reps > 0 and len(reps_with_mut) == total_reps else "; ".join(reps_with_mut),
                },
                seen_rows,
            )

    ref_seq = get_ref_seq(gene)
    if gene == "Cob":
        aa_ref_seq = nt_to_aa(ref_seq, mtDNA=True)
    else:
        aa_ref_seq = nt_to_aa(ref_seq)
    aligned_aa_seq, aligned_ref_aa_seq = align_seqs(aa_ref_seq, aa_seq)
    ref_alignment_map = alignment_mapping(aligned_ref_aa_seq)
    sample_alignment_map = alignment_mapping(aligned_aa_seq)

    # Check for novel mutations
    for i, (ref_aa, seq_aa) in enumerate(zip(aligned_ref_aa_seq, aligned_aa_seq), start=1):
        if ref_aa == seq_aa or seq_aa == "X":
            continue

        unaligned_ref_pos = ref_alignment_map[i - 1]
        unaligned_pos = sample_alignment_map[i - 1]

        is_known_variant = any(
            item["Gene"] == gene
            and unaligned_ref_pos == item["Posn"]
            and seq_aa in (item["Var1"], item["Var2"], item["Var3"])
            for item in variants
        )
        if is_known_variant:
            continue
        
        # Find which representatives have this novel mutation
        reps_with_mut = find_representatives_with_mutation(rep_sequences, gene, unaligned_ref_pos, seq_aa)
        total_reps = len(rep_sequences.get(gene, {}))
        
        add_row(
            matching_rows,
            {
                "Sample Name": str(snakemake.wildcards.sample),
                "Gene": gene,
                "WT amino acid": ref_aa,
                "MUT amino acid": seq_aa,
                "Ref. Organism": str(snakemake.wildcards.organism).title(),
                "Ref. aa position": unaligned_ref_pos,
                "Sample amino acid position": unaligned_pos,
                "Reference": "N/A",
                "Organism": str(snakemake.wildcards.organism).title(),
                "Representatives with mutation": "All representatives" if total_reps > 0 and len(reps_with_mut) == total_reps else "; ".join(reps_with_mut),
            },
            seen_rows,
        )

pd.DataFrame(matching_rows).to_csv(snakemake.output.multiqc_fung, index=None)

# Adding R-genes and Avr records to the multiqc report

def out_vir(geneID, sequence):
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
        new_row = out_vir(geneID, sequence)
        if new_row:
            avr_rows.append(new_row)
    elif geneID.startswith('Sr'):
        print(geneID)
        new_row = out_vir(geneID, sequence)
        if new_row:
            sr_rows.append(new_row)

pd.DataFrame(avr_rows).to_csv(snakemake.output.multiqc_avr, index=None)
pd.DataFrame(sr_rows).to_csv(snakemake.output.multiqc_rgene, index=None)