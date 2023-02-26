import pandas as pd

pd.concat([pd.read_csv(tsv, sep="\t") for tsv in snakemake.input.tsvs]).to_csv(
    snakemake.output.tsv, index=False, sep="\t"
)
