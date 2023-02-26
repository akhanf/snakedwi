import pandas as pd

df = pd.DataFrame()

missing_subj_zip_list = snakemake.params.missing_subject_zip_list

df["participant_id"] = missing_subj_zip_list["subject"]

if "session" in missing_subj_zip_list:
    df["session"] = missing_subj_zip_list["session"]

df.to_csv(snakemake.output.tsv, sep="\t", index=False)
