#!/usr/bin/env python

import pandas as pd


def create_missing_subj_tsv(
    missing_subj_list: dict,
    out_tsv: str,
) -> None:
    df = pd.DataFrame()
    df["participant_id"] = missing_subj_list["subject"]

    if "session" in missing_subj_list:
        df["session"] = missing_subj_list["session"]

    df.to_csv(out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    create_missing_subj_tsv(
        missing_subj_list=(snakemake.params.missing_subj_zip_list,),  # noqa: F821,
        out_tsv=snakemake.output.tsv,  # noqa: F821
    )
