import pandas as pd
from typing import Dict, IO, Iterable, List
import gzip
from Bio.SeqIO import parse


def file(path, mode='rt') -> IO:
    """Create a file object from path - infers compression from extension."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def read_pileup_depths(
    pileup_path: str,
    gene_lengths: Dict[str, int]
) -> Dict[str, List[int]]:

    depths = {
        gene: [0 for _ in range(length)]
        for gene, length in gene_lengths.items()
    }

    for line in file(pileup_path):
        gene, pos, _, depth, *_ = line.split()
        depths[gene][int(pos) - 1] = int(depth)

    return depths


def known_bases(seq: Iterable[str]) -> int:
    return sum(base != 'N' for base in seq)


def primer_performance(
    pileup_path: str,
    primers_path: str,
    consensus_path: str,
    primer_performance_out_path: str,
    thresholds: List[int] = [2, 10, 20],
):

    primers = pd.read_csv(primers_path)

    gene_lengths = primers.groupby('gene').gene_length.first().to_dict()
    depths = read_pileup_depths(pileup_path, gene_lengths)

    gene_thresholds = {}
    amplicon_thresholds = {}
    gene_n_used = {}
    amplicon_n_used = {}
    consensus = {r.id: str(r.seq) for r in parse(consensus_path, 'fasta')}
    
    for gene, gene_depths in depths.items():
        gene_seq = consensus[gene]
        gene_thresholds[gene] = {
            f'gene_ge_{threshold}': sum(d >= threshold for d in gene_depths)
            for threshold in thresholds
        }
        gene_n_used[gene] = [known_bases(gene_seq)]
        for _, amplicon in primers[primers.gene == gene].iterrows():
            amplicon_depths = gene_depths[amplicon.start:amplicon.end]
            amplicon_thresholds[amplicon.amplicon] = {
                f'amplicon_ge_{threshold}': sum(
                    d >= threshold for d in amplicon_depths
                )
                for threshold in thresholds
            }
            amplicon_seq = gene_seq[amplicon.start:amplicon.end]
            amplicon_n_used[amplicon.amplicon] = [known_bases(amplicon_seq)]

    primers = pd.merge(
        primers,
        pd.DataFrame(amplicon_thresholds).T.reset_index().rename(
            columns={'index': 'amplicon'}
        ),
        on='amplicon'
    )
    primers = pd.merge(
        primers,
        pd.DataFrame(gene_thresholds).T.reset_index().rename(
            columns={'index': 'gene'}
        ),
        on='gene'
    )
    primers = pd.merge(
        primers,
        pd.DataFrame(gene_n_used, index=['gene_n_used']).T.reset_index().rename(
            columns={'index': 'gene'}
        ),
        on='gene'
    )
    primers = pd.merge(
        primers,
        pd.DataFrame(amplicon_n_used, index=['amplicon_n_used']).T.reset_index().rename(
            columns={'index': 'amplicon'}
        ),
        on='amplicon'
    )

    for unit in ['amplicon', 'gene']:
        for threshold in thresholds:
            primers[f'pct_{unit}_ge_{threshold}'] = (
                100 * primers[f'{unit}_ge_{threshold}'] / primers[f'{unit}_length']
            ).round(4)
            primers[f'pct_{unit}_used'] = (
                100 * primers[f'{unit}_n_used'] / primers[f'{unit}_length']
            ).round(4)
            
    primers.to_csv(primer_performance_out_path, index=None)

if __name__ == '__main__':
    primer_performance(
        pileup_path=str(snakemake.input.pileup),
        primers_path=str(snakemake.input.primers),
        consensus_path=str(snakemake.input.consensus),
        primer_performance_out_path=str(snakemake.output.primer_performance),
        thresholds=snakemake.config['primer_performance_thresholds'],
    )