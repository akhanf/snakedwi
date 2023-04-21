#!/usr/bin/env python
import json

import numpy as np


def get_slspec_txt(in_jsons: str, out_txt: str) -> None:
    # get slice timing from first json file only.
    # TODO: Add check to ensure all are identical
    with open(in_jsons) as f:
        json_dwi = json.load(f)

    slice_timing = np.array(json_dwi["SliceTiming"])
    # Order in which slices are acquires
    arg_sorting = np.argsort(slice_timing)
    # Number of slice groups (by # of unique slice times)
    sg = len(set(slice_timing))
    # Multiband factor
    mb = int(len(slice_timing) / sg)
    # Reshape arg_sorting to get slspec and sort each row to make prettier
    slspec = np.sort(np.reshape(arg_sorting, [sg, mb]), axis=1)
    # Write to text
    np.savetxt(out_txt, slspec, fmt="%d")


if __name__ == "__main__":
    get_slspec_txt(
        in_json=snakemake.input.dwi_jsons,  # noqa: F821
        out_txt=snakemake.output.eddy_slspec_txt,  # noqa: F821
    )
