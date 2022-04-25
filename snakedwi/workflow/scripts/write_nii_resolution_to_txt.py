import numpy as np
import nibabel as nib

res_mm = nib.load(snakemake.input[0]).header.get_zooms()

with open(snakemake.output[0], "w") as f:
    f.write("x".join([str(vox) for vox in res_mm]) + "mm")
