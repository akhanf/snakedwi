#!/usr/bin/env python
import numpy as np


def get_synb0_phase_encode(in_fname: str, out_fname: str) -> None:
    # Read concatenated file
    with open(in_fname, "r", encoding="utf-8") as f:
        phenc_concat = np.loadtxt(f)

    # Update file for synb0 expected input
    # If multidimensional, just grab first row
    synb0_phenc = (
        np.array([phenc_concat[0], phenc_concat[0]])
        if phenc_concat.ndim > 1
        else np.array([phenc_concat, phenc_concat])
    )

    synb0_phenc[-1][-1] = 0

    # Write to output
    np.savetxt(out_fname, synb0_phenc)


if __name__ == "__main__":
    get_synb0_phase_encode(
        in_fname=snakemake.input.phenc_concat,  # noqa: F821
        out_fname=snakemake.output.synb0_phenc,  # noqa: F821
    )
