import json
import sys

for json_file in snakemake.input:

    with open(json_file, "r") as f:
        json_dwi = json.load(f)

        if (snakemake.config["slspec_txt"] == False) and (
            snakemake.config["eddy_no_s2v"] == False
        ):
            if "SliceTiming" not in json_dwi:
                print(f"ERROR: SliceTiming not found in {json_file}")
                print("You must do one of the following:")
                print(" 1. add the SliceTiming field to your dwi JSON files")
                print(" 2. use the --slspec_txt option to provide a global slspec file")
                print(
                    " 3. use the --eddy_no_s2v option to disable slice-to-volume correction"
                )

                sys.exit(1)

        if snakemake.config["no_topup"] == False:

            if "PhaseEncodingDirection" not in json_dwi:
                print(f"ERROR: PhaseEncodingDirection not found in {json_file}")
                print("You must do one of the following:")
                print(" 1. add the PhaseEncodingDirection field to your dwi JSON files")
                print(" 2. use the --no_topup option to disable top-up")
                sys.exit(1)

            if "EffectiveEchoSpacing" not in json_dwi:
                effechospc = snakemake.config["default_effective_echo_spacing"]
                print(f"WARNING: EffectiveEchoSpacing not found in {json_file}")
                print(f"If not defined, a default value of {effechospc} will be used")
