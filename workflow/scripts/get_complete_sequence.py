from Bio import SeqIO
import pysam

def read_index(ref_fnafai):
    contig_info = {}
    with open(ref_fnafai, 'r') as index_file:
        for line in index_file:
            fields = line.strip().split('\t')
            contig_name = fields[0]
            contig_length = int(fields[1])
            contig_info[contig_name] = contig_length
    return contig_info


def get_complete_sequence(bam_path, ref_path, index_path, output_path):
    ref_seqs = {record.id: list(record.seq) for record in SeqIO.parse(ref_path, "fasta")}
    contig_info = read_index(index_path)

    # Start with an empty sequence
    complete_sequences = {contig: ['N'] * length for contig, length in contig_info.items()}
    
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # Check reads in the bam file
    for read in bamfile.fetch():
        contig = read.reference_name
        if contig in complete_sequences:
            ref_seq = complete_sequences[contig]
            read_seq = read.query_sequence
            read_pos = read.reference_start

            # Update reference sequence with the read sequence
            # If the len(read) > len(reference), stop reading at the end of the reference
            for i, base in enumerate(read_seq):
                while read_pos + i < len(ref_seq):
                    ref_seq[read_pos + i] = base
                    break


    with open(output_path, "w") as output_file:
        for contig, seq in complete_sequences.items():
            output_file.write(f">{contig}\n{''.join(seq)}\n")

    bamfile.close()

get_complete_sequence(
    bam_path=snakemake.input.bam,
    ref_path=snakemake.input.ref,
    index_path=snakemake.input.index,
    output_path=snakemake.output.complete_seq,
)