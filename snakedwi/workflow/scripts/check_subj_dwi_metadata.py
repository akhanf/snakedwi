import json
import sys
import re
import pandas as pd
from snakemake.shell import shell

phase_encoding_directions = []
phase_encoding_axes = []

phenc_dirs = []

eddy_s2v = snakemake.config["use_eddy_s2v"]

df = pd.DataFrame()
row = dict()


has_slice_timing = True

for json_file in snakemake.input:

    with open(json_file, "r") as f:
        json_dwi = json.load(f)

    if eddy_s2v and (snakemake.config["slspec_txt"] == False):
        if "SliceTiming" not in json_dwi:
            has_slice_timing = False

            # print(f"WARNING: disabling s2v for {snakemake.wildcards}")

            eddy_s2v = False

            print(
                f"Eddy slice to volume enabled, but SliceTiming not found in {json_file}"
            )
            print("You must do one of the following:")
            print(" 1. add the SliceTiming field to your dwi JSON files")
            print(" 2. use the --slspec_txt option to provide a global slspec file")
            print(
                " 3. do not set the --use_eddy_s2v option to disable slice-to-volume correction"
            )

            sys.exit(1)

    if "PhaseEncodingDirection" not in json_dwi:

        if not snakemake.config["default_phase_encoding_direction"] == "":
            print(f"WARNING: setting default PhaseEncodingDirection")
            json_dwi["PhaseEncodingDirection"] = snakemake.config[
                "default_phase_encoding_direction"
            ]
        else:
            if "PhaseEncodingAxis" in json_dwi:
                print(
                    f"WARNING: assuming PhaseEncodingDirection from PhaseEncodingAxis"
                )
                json_dwi["PhaseEncodingDirection"] = json_dwi["PhaseEncodingAxis"]
            else:
                print(f"ERROR: PhaseEncodingDirection not found in {json_file}")
                print(
                    "You must add the PhaseEncodingDirection field to your dwi JSON files, or use the --default_phase_encoding_direction CLI option"
                )
                sys.exit(1)

    if "EffectiveEchoSpacing" in json_dwi:
        eff_echo=json_dwi["EffectiveEchoSpacing"]
    elif "EstimatedEffectiveEchoSpacing" in json_dwi:
        eff_echo=json_dwi["EstimatedEffectiveEchoSpacing"]
    else:
        print("EffectiveEchoSpacing not defined in JSON, using default value")
        eff_echo = snakemake.config[
            "default_effective_echo_spacing"
        ]


    phenc_dirs.append(json_dwi["PhaseEncodingDirection"])

    phase_encoding_axes.append(json_dwi["PhaseEncodingDirection"][0])
    if json_dwi["PhaseEncodingDirection"][-1] == "-":
        phase_encoding_directions.append("-")
    else:
        phase_encoding_directions.append("+")

# print(f"PhaseEncodingDirections: {phenc_dirs}")


if len(set(phase_encoding_axes)) > 1:
    #    print(f"WARNING: Multiple phase encoding axes used {phase_encoding_axes}")

    # select indices to use
    # pick LR/RL *only* if it also has multiple directions
    #  otherwise pick AP/PA

    # if we have LR+RL, use that axis:
    lr_indices = [i for i, pe in enumerate(phase_encoding_axes) if pe == "i"]
    ap_indices = [i for i, pe in enumerate(phase_encoding_axes) if pe == "j"]

    if len(set([phase_encoding_directions[i] for i in lr_indices])) == 2:
        # if we have lr+rl, use it
        use_indices = lr_indices
    else:
        # otherwise we use ap (even if we don't have multiple pedirs for ap)
        use_indices = ap_indices

else:
    # use all indices
    use_indices = [i for i in range(len(snakemake.input))]


# select the indices
phase_encoding_axes = [phase_encoding_axes[i] for i in use_indices]
phase_encoding_directions = [phase_encoding_directions[i] for i in use_indices]


# make the output folder
shell("mkdir -p {snakemake.output.workflowopts}")

# print(use_indices)

write_indices = ",".join([f"{ind}" for ind in use_indices])

# write the indices to file
shell("touch {snakemake.output.workflowopts}/indices-{write_indices}")

if len(set(phase_encoding_directions)) < 2:
    if phase_encoding_axes[0] == "i":

        print(
            "ERROR: Only one phase encoding direction with LR or RL - sdc will be inaccurate, so failing instead!"
        )
        sys.exit(1)

    if snakemake.config["use_syn_sdc"]:

        # print(
        #    f"Opposing phase encoding directions not available, {phase_encoding_directions}, using syn for sdc"
        # )
        shell("touch {snakemake.output.workflowopts}/sdc-syn")
    elif snakemake.config["use_synthsr_sdc"]:

        # print(
        #    f"Opposing phase encoding directions not available, {phase_encoding_directions}, using synthSR+Syn for sdc"
        # )
        shell("touch {snakemake.output.workflowopts}/sdc-synthsr")

    else:
        # print(
        #    f"Opposing phase encoding directions not available, {phase_encoding_directions}, skipping sdc"
        # )
        shell("touch {snakemake.output.workflowopts}/sdc-none")
else:
    # print(
    #    f"Opposing phase encoding directions are available, {phase_encoding_directions}, using topup for sdc"
    # )
    shell("touch {snakemake.output.workflowopts}/sdc-topup")

if eddy_s2v:
    # print("Enabling eddy s2v in the workflow")
    shell("touch {snakemake.output.workflowopts}/eddys2v-yes")
else:
    # print("Disabling eddy s2v in the workflow")
    shell("touch {snakemake.output.workflowopts}/eddys2v-no")

# print("Writing phase encoding axis")
shell("touch {snakemake.output.workflowopts}/PEaxis-{phase_encoding_axes[0]}")

# sets the index column as sub-{subject} or sub-{subject}_ses-{session}
row[snakemake.params.index_col_name] = [snakemake.params.index_col_value]
row["has_slice_timing"] = [{True: "yes", False: "no"}[has_slice_timing]]
row["pedirs"] = [
    ",".join(
        [
            f"{peax}{pedir}"
            for peax, pedir in zip(phase_encoding_axes, phase_encoding_directions)
        ]
    )
]

df = pd.DataFrame.from_dict(row)
# write to output file
df.to_csv(snakemake.output.metadata, sep="\t", index=False)
