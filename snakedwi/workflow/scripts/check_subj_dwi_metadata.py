import json
import sys
import re
from snakemake.shell import shell

phase_encoding_directions = []
phase_encoding_axes = []

phenc_dirs = []

for json_file in snakemake.input:

    with open(json_file, "r") as f:
        json_dwi = json.load(f)

    if (snakemake.config["slspec_txt"] == False) and (
        snakemake.config["use_eddy_s2v"] == True
    ):
        if "SliceTiming" not in json_dwi:
            print(
                f"ERROR: Eddy slice to volume enabled, but SliceTiming not found in {json_file}"
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

if len(set(phase_encoding_directions)) < 2:
    print(
        f"Opposing phase encoding directions not available, {phase_encoding_directions}, using syn for sdc"
    )
    shell("mkdir -p {snakemake.output} && touch {snakemake.output}/sdc-syn")
else:
    print(
        f"Opposing phase encoding directions are available, {phase_encoding_directions}, using topup for sdc"
    )
    shell("mkdir -p {snakemake.output} && touch {snakemake.output}/sdc-topup")
