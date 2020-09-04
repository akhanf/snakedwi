import json
import nibabel as nib
import numpy as np


with open(snakemake.input.shells) as f:
  shells_dict = json.load(f)

#get bval parameter:
bval = snakemake.params.bval

#input dwi
dwi_nib = nib.load(snakemake.input.dwi)
print(dwi_nib.shape)

#create output shape
newshape = np.array(dwi_nib.shape[:3])

avg_shell = np.zeros(newshape.astype(int))


indices = shells_dict['shell_to_vol'][bval] 
avg_shell = np.mean(dwi_nib.get_fdata()[:,:,:,indices],3)

#now save as image
avg_shell_nii = nib.Nifti1Image(avg_shell, affine=dwi_nib.affine )
avg_shell_nii.to_filename(snakemake.output[0])
