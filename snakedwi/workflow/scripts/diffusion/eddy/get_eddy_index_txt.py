#!/usr/bin/env python
from typing import List

import nibabel as nib
import numpy as np


def get_eddy_index_txt(in_niis: List[str], out_txt: str) -> None:
    # Get each nifti's imsize
    imsizes = [nib.load(nii).header.get_data_shape() for nii in in_niis]

    eddy_idxes = []
    for i in range(len(imsizes)):
        index = i + 1
        if len(imsizes[i]) < 4:  # 3D vol
            eddy_idxes.append(index)
        else:
            for j in range(imsizes[i][3]):
                eddy_idxes.append(index)

    # Save to file
    np.savetxt(out_txt, np.array(eddy_idxes), fmt="%d")


if __name__ == "__main__":
    get_eddy_index_txt(
        in_niis=snakemake.input.dwi_niis,  # noqa: F821
        out_txt=snakemake.output.eddy_index_txt,  # noqa: F821
    )
