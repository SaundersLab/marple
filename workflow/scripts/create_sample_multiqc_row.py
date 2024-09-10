import pandas as pd
from Bio.SeqIO import parse
import os

tree_input_path = str(snakemake.input.tree_input)
tree_input_seq = next(str(r.seq) for r in parse(tree_input_path, 'fasta'))
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