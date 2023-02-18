import nibabel as nib
import json
from sdcflows.transform import disp_to_fmap
import numpy as np


# adapted from sdcflows DisplacementsField2Fieldmap Interface
with open(snakemake.input.dwi_json, "r") as f:
    json_dwi = json.load(f)

# get phenc dir, and calculate required readout time
pe_dir = json_dwi["PhaseEncodingDirection"]
shape_dwi = nib.load(snakemake.input.dwi_nii).get_fdata().shape
n_pe = shape_dwi["ijk".index(pe_dir[0])]
ro_time = json_dwi["EffectiveEchoSpacing"] * (n_pe - 1)

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
