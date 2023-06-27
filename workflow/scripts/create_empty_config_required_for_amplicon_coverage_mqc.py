with open(str(snakemake.output.placeholder_mqc), 'w') as f:
    f.write(f'id: "{snakemake.wildcards.sample}_placeholder"\n')
