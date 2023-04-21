#!/usr/bin/env python
from typing import List

import pandas as pd


def concat_subj_metadata(in_tsvs: List[str], out_tsv: str) -> None:
    pd.concat([pd.read_csv(tsv, sep="\t") for tsv in in_tsvs]).to_csv(
        out_tsv, sep="\t", index=False
    )


if __name__ == "__main__":
    concat_subj_metadata(
        in_tsvs=snakemake.input.tsvs,  # noqa: F821
        out_tsv=snakemake.output.tsv,  # noqa: F821
    )
