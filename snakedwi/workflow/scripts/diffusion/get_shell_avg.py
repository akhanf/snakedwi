#!/usr/bin/env python
import json
import sys

import nibabel as nib
import numpy as np


def get_avg_shell(dwi: str, shells: str, bval: str, out_nii: str) -> None:
    # Load relevant data
    dwi_nii = nib.load(dwi)
    dwi_data, dwi_affine = dwi_nii.get_fdata(), dwi_nii.affine
    dwi_shape = dwi_nii.shape

    with open(shells) as shell_fname:
        shells_dict = json.load(shell_fname)

    # Create output shape
    newshape = np.array(dwi_shape[:3])
    avg_shell = np.zeros(newshape.astype(int))

    indices = shells_dict["shell_to_vol"][bval]

    # 3D vol (e.g. b0 only) - grab directly
    if len(dwi_shape) == 3 and len(indices) == 1 and indices[0] == 0:
        avg_shell = dwi_data[:, :, :]
    # 4D vol - extract indices and average
    elif len(dwi_shape) == 4:
        avg_shell = np.mean(dwi_data[:, :, :, indices], 3)
    # Something weird with indices and volume
    else:
        print("Unable to get map indices to get avg shell...exiting!")
        sys.exit(1)

    # Save output
    avg_shell_nii = nib.Nifti1Image(avg_shell, affine=dwi_affine)
    avg_shell_nii.to_filename(out_nii)


if __name__ == "__main__":
    get_avg_shell(
        dwi=snakemake.input.dwi,  # noqa: F821
        shells=snakemake.input.shells,  # noqa: F821
        bval=snakemake.params.bval,  # noqa: F821
        out_nii=snakemake.output.avgshell,  # noqa: F821
    )
