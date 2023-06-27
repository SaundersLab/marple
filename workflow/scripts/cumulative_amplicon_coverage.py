import pandas as pd

primers_performance = pd.read_csv(str(snakemake.input.primer_performance))

with open(str(snakemake.output.cumulative_amplicon_coverage), 'w') as f:
    for pct in range(101):
        n = primers_performance[primers_performance.amplicon_pct_used >= pct].shape[0]
        f.write(f'{pct},{n}\n')
