import json
import sys

import nibabel as nib
import numpy as np
from sdcflows.transform import disp_to_fmap


def displacement_field_to_fmap(
    disp_field: str,
    dwi_nii: str,
    dwi_json: str,
    demean: bool,
    out_fmap: str,
    smk_config: dict,
) -> None:
    """Adapted from SDCFlows DisplacementField2Fieldmap Interface"""

    with open(dwi_json, "r") as f:
        json_dwi = json.load(f)

    # Get phenc dir & calculate required readout time
    if "PhaseEncodingDirection" in json_dwi:
        pe_dir = json_dwi["PhaseEncodingDirection"]
    else:
        if not smk_config["default_phase_encoding_direction"] == "":
            print("WARNING: setting default PhaseEncodingDirection")
            pe_dir = smk_config["default_phase_encoding_direction"]
        else:
            if "PhaseEncodingAxis" in json_dwi:
                print(
                    "WARNING: assuming PhaseEncodingDirection from "
                    "PhaseEncodingAxis"
                )
                pe_dir = json_dwi["PhaseEncodingAxis"]
            else:
                print(
                    "ERROR: PhaseEncodingDirection not found in dwi JSON.\n"
                    "You must add the PhaseEncodingDirection field to your "
                    "dwi JSON files, or use the "
                    "--default_phase_encoding_direction CLI option"
                )
                sys.exit(1)

    if "EffectiveEchoSpacing" in json_dwi:
        eff_echo = json_dwi["EffectiveEchoSpacing"]
    elif "EstimatedEffectiveEchoSpacing" in json_dwi:
        eff_echo = json_dwi["EstimatedEffectiveEchoSpacing"]
    else:
        print("EffectiveEchoSpacing not defined in JSON, using default value")
        eff_echo = smk_config["default_effective_echo_spacing"]

    dwi_shape = nib.load(dwi_nii).get_fdata().shape
    n_pe = dwi_shape["ijk".index(pe_dir[0])]
    ro_time = eff_echo * (n_pe - 1)

    # pass this on to the sdcflows function
    nib_fmap = disp_to_fmap(nib.load(disp_field), ro_time, pe_dir)

    # demean the fieldmap
    if demean:
        data = np.asanyarray(nib_fmap.dataobj)
        data -= np.median(data)

    _ = nib_fmap.__class__(
        data.astype("float32"),
        nib_fmap.affine,
        nib_fmap.header,
    )

    # Save fmap
    nib_fmap.to_filename(out_fmap)


if __name__ == "__main__":
    displacement_field_to_fmap(
        dwi_nii=snakemake.input.dwi_nii,  # noqa: F821
        dwi_json=snakemake.input.dwi_json,  # noqa: F821
        disp_field=snakemake.input.disp_field,  # noqa: F821
        demean=snakemake.params.demean,  # noqa: F821
        out_fmap=snakemake.output.fmap,  # noqa: F821
        smk_config=snakemake.config,  # noqa: F821
    )
