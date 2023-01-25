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
    phenc_axis = json_dwi["PhaseEncodingDirection"][0]
    phenc_string = "PhaseEncodingDirection"

# NOTE: PhaseEncodingAxis doesn't encode + or -, so if this is used then you cannot use rev phase encoding distortion correction (but it can still be used if there is just a single scan)
elif "PhaseEncodingAxis" in json_dwi:
    phenc_axis = json_dwi["PhaseEncodingAxis"][0]  # either i, j, k
    phenc_string = "PhaseEncodingAxis"

if phenc_axis == "i":
    vec = np.array([1, 0, 0])
elif phenc_axis == "j":
    vec = np.array([0, 1, 0])
elif phenc_axis == "k":
    vec = np.array([0, 0, 1])

# print(f'vec: {vec}')
# print(f'imsize: {imsize}')f

numPhaseEncodes = imsize[np.where(vec > 0)]

# print(f'numPhaseEncodes: {numPhaseEncodes}')

# check for i-, j-, k-; flip to -1 if so..
if len(json_dwi[phenc_string]) == 2:
    if json_dwi[phenc_string][1] == "-":
        vec[np.where(vec > 0)] = -1

if "EffectiveEchoSpacing" in json_dwi:
    phenc_line = np.hstack(
        [vec, np.array(json_dwi["EffectiveEchoSpacing"] * numPhaseEncodes)]
    )
elif "EstimatedEffectiveEchoSpacing" in json_dwi:
    print("Estimtated Effective Echo Spacing found, using that")
    phenc_line = np.hstack(
        [vec, np.array(json_dwi["EstimatedEffectiveEchoSpacing"] * numPhaseEncodes)]
    )
else:
    print("EffectiveEchoSpacing not defined in JSON, using default value")
    json_dwi["EffectiveEchoSpacing"] = snakemake.config[
        "default_effective_echo_spacing"
    ]
    phenc_line = np.hstack(
        [vec, np.array(json_dwi["EffectiveEchoSpacing"] * numPhaseEncodes)]
    )


# create the phenc_line row
# phenc_line = np.hstack(
#    [vec, np.array(json_dwi["EffectiveEchoSpacing"] * numPhaseEncodes)]
# )

# replicate to the number of volumes, if it is 4d
if len(imsize) == 4:
    phenc_data = np.tile(phenc_line, [imsize[3], 1])
else:
    phenc_data = np.column_stack(phenc_line)
    # need to column_stack so that it becomes 2d array
    # otherwwise savetxt will always save 1d array as column..


# save to txt
np.savetxt(snakemake.output.phenc_txt, phenc_data, fmt="%.5f")
