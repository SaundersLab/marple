from re import sub, search
from typing import Tuple, Dict, Union, List, IO
from Bio.SeqIO import parse
import gzip

def file(path, mode='rt') -> IO:
    """Create a file object from path - infers compression from extension."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def base_ratios_from_reads(ref: str,
    depth: int,
    reads: str,
) -> Dict[str, float]:

    # Don't count mapping qualities at the start of read segments as bases
    # e.g. ^G. is a read segment start and the ASCII minus 33 of G is the
    # mapping quality. Drop the 2 characters that are not reads.
    if '^' in reads:
        reads = sub('\^.', '', reads)

    # Remove pileup indels so you don't count the bases in an indel e.g. +3ATT
    if '+' in reads or '-' in reads:
        while match := search(r'[+-](\d+)', reads):
            reads = reads[:match.start()] + reads[match.end() + int(match[1]):]

    upper_reads = reads.upper()
    base_counts = {base: upper_reads.count(base) for base in 'ACTG'}

    upper_ref = ref.upper()
    if upper_ref in 'ACTG':
        base_counts[upper_ref] += reads.count('.') + reads.count(',')

    # Only calculate ratios for present bases
    return {base: n / depth for base, n in base_counts.items() if n}


def get_genotype_and_valid_base_ratios(
    ref: str,
    base_ratios: Dict[str, float],
    hetero_min: float,
    hetero_max: float,
) -> Union[Tuple[str, Dict[str, float]], Tuple[None, None]]:
    """Use the heterozygosity thresholds to filter the base ratios and determine the genotype:
    Genotype	Condition	                                        Example
    0/0	        pref ≥ hetero_max or only refbase qualified	        G GG 0/0 1
    0/1	        ref base and 1 alt base qualified	                T CT 0/1 0.217,0.783
    1/1	        palt 1 ≥ hetero_max or only 1 alt base qualified	T CC 1/1 0.822
    1/2	        2 alt bases qualified and ref base did not	        C AT 1/2 0.476,0.524
    ?	        3+ bases qualified	                                T ACT ? 0.2,0.6,0.2
    None        No bases qualify
    """
    assert hetero_max > hetero_min

    valid_base_ratios = {}
    for base, ratio in base_ratios.items():
        if ratio >= hetero_min:
            valid_base_ratios[base] = ratio
            if ratio >= hetero_max:
                return '0/0' if base == ref else '1/1', valid_base_ratios

    if not valid_base_ratios:
        return None, None

    # If only ref has min <= freq <= max
    if list(valid_base_ratios) == [ref]:
        return None, None

    if len(valid_base_ratios) == 1:
        genotype = '1/1'
    elif len(valid_base_ratios) == 2:
        genotype = '0/1' if ref in valid_base_ratios else '1/2'
    else:
        genotype = '?'

    return genotype, valid_base_ratios


def pileup_row_to_consensus(
    row: str,
    snp_file_handle: IO,
    min_ref_depth: int = 2,
    min_snp_depth: int = 10,
    hetero_min: float = .2,
    hetero_max: float = .8,
) -> Tuple[str, int, str]:
    alleles_to_code = {
        'A': 'A',  'C': 'C',  'G': 'G',  'T': 'T',
        'AT': 'W',  'CG': 'S',  'AC': 'M',  'GT': 'K',  'AG': 'R',  'CT': 'Y',
        'TA': 'W',  'GC': 'S',  'CA': 'M',  'TG': 'K',  'GA': 'R',  'TC': 'Y',
    }
    contig, pos_str, ref, depth_str, reads, *_ = row.split('\t')
    depth = int(depth_str)
    pos = int(pos_str)
    null_consensus = contig, pos - 1, 'N'
    if depth < min_ref_depth:
        return null_consensus
    base_ratios = base_ratios_from_reads(ref, depth, reads)
    if list(base_ratios) == [ref]:
        return contig, pos - 1, ref
    if depth < min_snp_depth:
        return null_consensus
    genotype, valid_base_ratios = get_genotype_and_valid_base_ratios(
        ref, base_ratios, hetero_min, hetero_max
    )
    if valid_base_ratios is None:
        return null_consensus
    
    alleles = ''.join(valid_base_ratios)
    
    # write the SNP to a file
    if snp_file_handle is not None:
        ratios = ','.join(map(lambda v: str(round(v, 3)), base_ratios.values()))
        snp_row = (contig, pos - 1, ref, depth, alleles, ratios, genotype)
        snp_file_handle.write('\t'.join(map(str, snp_row)) + '\n')

    if genotype == '?':
        return null_consensus
    
    return contig, pos - 1, alleles_to_code[alleles]


def pileup_to_consensus(
    pileup_path: str,
    ref_path: str,
    out_path: str = None,
    snp_file_path: str = None,
    **kwargs,
) -> Dict[str, List[str]]:

    # Start with an empty consensus
    consensus = {r.id: ['N'] * len(r.seq) for r in parse(ref_path, 'fasta')}
    with file(pileup_path) as f_in, file(snp_file_path, 'wt') as f_out:
        f_out.write('seqid\tpos\tref\tdepth\talleles\tratios\tgenotype\n')
        for row in f_in:
            contig, pos, consensus[contig][pos] = pileup_row_to_consensus(
                row=row, snp_file_handle=f_out, **kwargs
            )

    with file(out_path, 'wt') as f:
        for contig in sorted(consensus):
            codes = consensus[contig]
            f.write('>' + contig + '\n' + ''.join(codes) + '\n')
    return consensus

if __name__ == '__main__':
    pileup_to_consensus(
        pileup_path=str(snakemake.input.pileup),
        ref_path=str(snakemake.input.ref),
        out_path=str(snakemake.output.consensus),
        snp_file_path=str(snakemake.output.snps),
        min_ref_depth=int(snakemake.config['min_ref_depth']),
        min_snp_depth=int(snakemake.config['min_snp_depth']),
        hetero_min=float(snakemake.config['hetero_min']),
        hetero_max=float(snakemake.config['hetero_max']),
    )
