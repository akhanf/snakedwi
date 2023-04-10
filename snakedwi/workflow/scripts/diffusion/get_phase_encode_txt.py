#!/usr/bin/env python
import json
import sys

import nibabel as nib
import numpy as np


def get_phase_encode_txt(
    bzero_nii: str,
    json_file: str,
    phenc_txt: str,
    smk_config: dict,
) -> None:
    # Load relevant data
    with open(json_file) as json_fname:
        json_dwi = json.load(json_fname)
    bzero = nib.load(bzero_nii)
    imsize = np.array(bzero.header.get_data_shape())

    # Check metadata
    if "PhaseEncodingDirection" in json_dwi:
        phenc_dir = json_dwi["PhaseEncodingDirection"]
    else:
        if not smk_config["default_phase_encoding_direction"] == "":
            print("WARNING: setting default PhaseEncodingDirection")
            phenc_dir = smk_config["default_phase_encoding_direction"]
        else:
            if "PhaseEncodingAxis" in json_dwi:
                print(
                    "WARNING: assuming PhaseEncodingDirection from "
                    "PhaseEncodingAxis"
                )
                phenc_dir = json_dwi["PhaseEncodingAxis"]
            else:
                print(
                    "ERROR: PhaseEncodingDirection not found in dwi JSON.\n"
                    "You must add the PhaseEncodingDirection to your dwi JSON "
                    "files, or use the --default_phase_encoding_direction CLI "
                    "option."
                )
                sys.exit(1)

    if "EffectiveEchoSpacing" in json_dwi:
        eff_echo = json_dwi["EffectiveEchoSpacing"]
    elif "EstimatedEffectiveEchoSpacing" in json_dwi:
        eff_echo = json_dwi["EstimatedEffectiveEchoSpacing"]
    else:
        print("EffectiveEchoSpacing not defined in JSON, using default value")
        eff_echo = smk_config["default_effective_echo_spacing"]

    # Determine PE direction
    vecs = {
        "i": np.array([1, 0, 0]),
        "j": np.array([0, 1, 0]),
        "k": np.array([0, 0, 1]),
    }
    vec = vecs[phenc_dir[0]]
    # Compute number of phase encodes
    numPhaseEncodes = imsize[np.where(vec > 0)]
    # Check for i-, j-, k-; flip to -1 if so.
    if len(phenc_dir) == 2 and phenc_dir[1] == "-":
        vec[np.where(vec > 0)] = -1

    # Line to be added to out_file
    phenc_line = np.hstack([vec, np.array(eff_echo * numPhaseEncodes)])

    # If 4D - replicate to the number of volumes, if it is 4d
    # Otherwise column stack (2D) for savetxt (will always save 1D as column)
    if len(imsize) == 4:
        phenc_data = np.tile(phenc_line, [imsize[3], 1])
    else:
        phenc_data = np.column_stack(phenc_line)

    # Save to output file
    np.savetxt(phenc_txt, phenc_data, fmt="%.5f")


if __name__ == "__main__":
    get_phase_encode_txt(
        bzero_nii=snakemake.input.bzero_nii,  # noqa: F821
        json_file=snakemake.input.json,  # noqa: F821
        phenc_txt=snakemake.output.phenc_txt,  # noqa: F821
        smk_config=snakemake.config,  # noqa: F821
    )
