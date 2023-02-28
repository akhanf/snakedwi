import nibabel as nib
import json
from sdcflows.transform import disp_to_fmap
import numpy as np


# adapted from sdcflows DisplacementsField2Fieldmap Interface
with open(snakemake.input.dwi_json, "r") as f:
    json_dwi = json.load(f)

# get phenc dir, and calculate required readout time
if "PhaseEncodingDirection" in json_dwi:
    pe_dir = json_dwi["PhaseEncodingDirection"]
else:
    if not snakemake.config["default_phase_encoding_direction"] == "":
        print(f"WARNING: setting default PhaseEncodingDirection")
        pe_dir = snakemake.config[
            "default_phase_encoding_direction"
        ]
    else:
        if "PhaseEncodingAxis" in json_dwi:
            print(
                f"WARNING: assuming PhaseEncodingDirection from PhaseEncodingAxis"
            )
            pe_dir = json_dwi["PhaseEncodingAxis"]
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



shape_dwi = nib.load(snakemake.input.dwi_nii).get_fdata().shape
n_pe = shape_dwi["ijk".index(pe_dir[0])]
ro_time = ecff_echo * (n_pe - 1)

# pass this on to the sdcflows function
nib_fmap = disp_to_fmap(nib.load(snakemake.input.disp_field), ro_time, pe_dir)

# demean the fieldmap
if snakemake.params.demean:
    data = np.asanyarray(nib_fmap.dataobj)
    data -= np.median(data)

    fmapnii = nib_fmap.__class__(
        data.astype("float32"),
        nib_fmap.affine,
        nib_fmap.header,
    )

# save it
nib_fmap.to_filename(snakemake.output.fmap)
