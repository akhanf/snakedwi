#!/usr/bin/env python
import json
import sys

import nibabel as nib
import numpy as np


def get_shell_vols(dwi: str, shells: str, bval: str, out_nii: str) -> None:
    # Load relevant data
    dwi_nii = nib.load(dwi)
    dwi_data, dwi_affine = dwi_nii.get_fdata(), dwi_nii.affine
    dwi_shape = dwi_nii.shape

    with open(shells) as shell_fname:
        shells_dict = json.load(shell_fname)
    indices = shells_dict["shell_to_vol"][bval]

    # Create output shape, with 4th dim = # of volumes at specified shell
    newshape = np.array(
        [dwi_shape[0], dwi_shape[1], dwi_shape[2], len(indices)]
    )
    shell_vols = np.zeros(newshape.astype(int))

    # 3D vol (e.g. b0 only) - grab directly
    if len(dwi_shape) == 3 and len(indices) == 1 and indices[0] == 0:
        shell_vols = dwi_data[:, :, :]
    # 4D vol - extract indices
    elif len(dwi_shape) == 4:
        shell_vols = dwi_data[:, :, :, indices]
    # Something weird with indices and volume
    else:
        print("Unable to get map indices to get avg shell...exiting!")
        sys.exit(1)

    # Save output
    shell_vols_nii = nib.Nifti1Image(shell_vols, affine=dwi_affine)
    shell_vols_nii.to_filename(out_nii)


if __name__ == "__main__":
    get_shell_vols(
        dwi=snakemake.input.dwi,  # noqa: F821
        shells=snakemake.input.shells,  # noqa: F821
        bval=snakemake.params.bval,  # noqa: F821
        out_nii=snakemake.output.shell_vols,  # noqa: F821
    )
