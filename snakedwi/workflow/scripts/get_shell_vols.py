import json
import nibabel as nib
import numpy as np
import sys

with open(snakemake.input.shells) as f:
    shells_dict = json.load(f)

# get bval parameter:
bval = snakemake.params.bval

indices = shells_dict["shell_to_vol"][bval]

# input dwi
dwi_nib = nib.load(snakemake.input.dwi)

# create output shape, 4th dim length of indices (# of volumes at this shell)
newshape = np.array(
    [dwi_nib.shape[0], dwi_nib.shape[1], dwi_nib.shape[2], len(indices)]
)

shell_vols = np.zeros(newshape.astype(int))


if len(dwi_nib.shape) == 3 and len(indices) == 1 and indices[0] == 0:
    # we have 3d vol (e.g. b0 only), so just grab it..
    shell_vols = dwi_nib.get_fdata()[:, :, :]
elif len(dwi_nib.shape) == 4:
    # otherwise, pick out indices and average
    shell_vols = dwi_nib.get_fdata()[:, :, :, indices]
else:
    # if not either of these cases, then something weird with indices and volumes
    print("unable to get map indices to get shell volumes")
    sys.exit()

# now save as image
shell_vols_nii = nib.Nifti1Image(shell_vols, affine=dwi_nib.affine)
shell_vols_nii.to_filename(snakemake.output[0])
