#!/usr/bin/env python
import nibabel as nib


def write_nii_resolution(in_nii: str, out_txt: str) -> None:
    res_mm = nib.load(in_nii).header.get_zooms()

    with open(out_txt, "w") as f:
        f.write("x".join([str(vox) for vox in res_mm]) + "mm")


if __name__ == "__main__":
    write_nii_resolution(
        in_nii=snakemake.input.nii,  # noqa: F821
        out_txt=snakemake.output.txt_fname,  # noqa: F821
    )
