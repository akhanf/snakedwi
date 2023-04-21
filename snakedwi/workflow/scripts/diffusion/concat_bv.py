#!/usr/bin/env python
from typing import List

import numpy as np


def concat_bv(bv_files: List[str], out_bv: str) -> None:
    stacked_bv = np.hstack(
        [np.loadtxt(bv_file, ndmin=2) for bv_file in bv_files]
    )
    np.savetxt(out_bv, stacked_bv, fmt="%.5f")


if __name__ == "__main__":
    concat_bv(
        bv_files=snakemake.input.bv_files,  # noqa: F821
        out_bv=snakemake.output.out_fname,  # noqa: F821
    )
