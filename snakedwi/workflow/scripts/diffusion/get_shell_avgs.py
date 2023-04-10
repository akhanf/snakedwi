#!/usr/bin/env python
import json

import nibabel as nib
import numpy as np


def get_shell_avgs(dwi: str, shells: str, out_nii: str) -> None:
    # Load relevant data
    dwi_nii = nib.load(dwi)
    dwi_data, dwi_affine = dwi_nii.get_fdata(), dwi_nii.affine

    with open(shells) as shell_fname:
        shells_dict = json.load(shell_fname)

    # Create output shape
    newshape = np.zeros(
        [
            4,
        ]
    )
    newshape[:3] = np.array(dwi_nii.shape[:3])
    newshape[3] = len(shells_dict["shells"])
    avg_shells = np.zeros(newshape.astype(int))

    # Extract and compute shell averages
    for idx, shell in enumerate(shells_dict["shells"]):
        indices = shells_dict["shell_to_vol"][str(shell)]
        avg_shells[:, :, :, idx] = np.mean(dwi_data[:, :, :, indices], 3)

    # Save output
    avg_shells_nii = nib.Nifti1Image(avg_shells, affine=dwi_affine)
    avg_shells_nii.to_filename(out_nii)


if __name__ == "__main__":
    get_shell_avgs(
        dwi=snakemake.input.dwi,  # noqa: F821
        shells=snakemake.input.shells,  # noqa: F821
        out_nii=snakemake.output.avgshells,  # noqa: F821
    )
