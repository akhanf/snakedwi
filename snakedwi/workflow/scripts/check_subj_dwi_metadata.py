import json
import sys
import re
from snakemake.shell import shell

phase_encoding_directions = []
phase_encoding_axes = []

phenc_dirs = []

eddy_s2v = snakemake.config["use_eddy_s2v"]


for json_file in snakemake.input:

    with open(json_file, "r") as f:
        json_dwi = json.load(f)

    if eddy_s2v and (snakemake.config["slspec_txt"] == False):
        if "SliceTiming" not in json_dwi:
            print(f"WARNING: disabling s2v for {snakemake.wildcards}")

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
        print(f"ERROR: PhaseEncodingDirection not found in {json_file}")
        print("You must add the PhaseEncodingDirection field to your dwi JSON files")
        sys.exit(1)

    phenc_dirs.append(json_dwi["PhaseEncodingDirection"])

    phase_encoding_axes.append(json_dwi["PhaseEncodingDirection"][0])
    if json_dwi["PhaseEncodingDirection"][-1] == "-":
        phase_encoding_directions.append("-")
    else:
        phase_encoding_directions.append("+")

print(f"PhaseEncodingDirections: {phenc_dirs}")

if len(set(phase_encoding_axes)) > 1:
    print(f"ERROR: Multiple phase encoding axes used {phase_encoding_axes}")
else:
    print(f"Phase encoding axes: {phase_encoding_axes}")


# make the output folder
shell("mkdir -p {snakemake.output}")

if len(set(phase_encoding_directions)) < 2:
    if snakemake.config["use_syn_sdc"]:

        print(
            f"Opposing phase encoding directions not available, {phase_encoding_directions}, using syn for sdc"
        )
        shell("touch {snakemake.output}/sdc-syn")
    else:
        print(
            f"Opposing phase encoding directions not available, {phase_encoding_directions}, skipping sdc"
        )
        shell("touch {snakemake.output}/sdc-none")
else:
    print(
        f"Opposing phase encoding directions are available, {phase_encoding_directions}, using topup for sdc"
    )
    shell("touch {snakemake.output}/sdc-topup")

if eddy_s2v:
    print("Enabling eddy s2v in the workflow")
    shell("touch {snakemake.output}/eddys2v-yes")
else:
    print("Disabling eddy s2v in the workflow")
    shell("touch {snakemake.output}/eddys2v-no")
