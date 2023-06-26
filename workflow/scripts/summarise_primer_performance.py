import pandas as pd
import os
from typing import List

def summarise_primer_performance(
    primer_performance_path: str,
    groupby_col: str,
    primer_performance_summary_path: List[str],
):
    primers = pd.read_csv(primer_performance_path)
    units = ['amplicon', 'gene']
    
    groupby = primers.groupby(groupby_col)
    unit_summaries = {}
    for unit in units:
        unit_summary = pd.DataFrame({
            col: groupby[col].sum()
            for col in primers
            if col.startswith(f'{unit}_') and ('pct' not in col)
        })
        for col in unit_summary:
            if 'length' not in col:
                unit_summary[f'{col.replace("_n_", "_pct_")}'] = (
                    100 * unit_summary[col] / unit_summary[f'{unit}_length']
                ).round(4)
        unit_summaries[unit] = unit_summary
    group_summary = pd.concat([unit_summaries[unit] for unit in units], axis=1)
    group_summary.reset_index().to_csv(primer_performance_summary_path, index=None)

if __name__ == '__main__':
    summarise_primer_performance(
        primer_performance_path=str(snakemake.input.primer_performance),
        groupby_col=snakemake.wildcards.groupby,
        primer_performance_summary_path=str(snakemake.output.primer_performance_summary),
    )