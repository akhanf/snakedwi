import nibabel as nib
import json
import numpy as np

# load nifti
bzero = nib.load(snakemake.input.bzero_nii)

# load json
with open(snakemake.input.json) as f:
    json_dwi = json.load(f)

# print(bzero.header)
# print(json_dwi)

imsize = np.array(bzero.header.get_data_shape())

if "PhaseEncodingDirection" in json_dwi:
    phenc_dir = json_dwi["PhaseEncodingDirection"]
else:
    if not snakemake.config["default_phase_encoding_direction"] == "":
        print(f"WARNING: setting default PhaseEncodingDirection")
        phenc_dir = snakemake.config[
            "default_phase_encoding_direction"
        ]
    else:
        if "PhaseEncodingAxis" in json_dwi:
            print(
                f"WARNING: assuming PhaseEncodingDirection from PhaseEncodingAxis"
            )
            phenc_dir = json_dwi["PhaseEncodingAxis"]
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



phenc_axis = phenc_dir[0]


if phenc_axis == "i":
    vec = np.array([1, 0, 0])
elif phenc_axis == "j":
    vec = np.array([0, 1, 0])
elif phenc_axis == "k":
    vec = np.array([0, 0, 1])


numPhaseEncodes = imsize[np.where(vec > 0)]


# check for i-, j-, k-; flip to -1 if so..
if len(phenc_dir) == 2:
    if phenc_dir[1] == "-":
        vec[np.where(vec > 0)] = -1


phenc_line = np.hstack(
    [vec, np.array(eff_echo * numPhaseEncodes)]
)


# replicate to the number of volumes, if it is 4d
if len(imsize) == 4:
    phenc_data = np.tile(phenc_line, [imsize[3], 1])
else:
    phenc_data = np.column_stack(phenc_line)
    # need to column_stack so that it becomes 2d array
    # otherwwise savetxt will always save 1d array as column..


# save to txt
np.savetxt(snakemake.output.phenc_txt, phenc_data, fmt="%.5f")
