#!/usr/bin/env python
import json

import numpy as np
from numpy import ndarray


def find_shells(bvals: ndarray):
    """Use histogram method to find number of shells, with anywhere from 10 to
    100 bins. Want to use the highest number of bins that doesn't split up the
    shell, such that the bin centers can be used as the shell bvalues."""
    shells = []
    for idx, nbins in enumerate(np.arange(10, 100, 5)):
        counts, bin_edges = np.histogram(bvals, bins=nbins)
        bin_lhs, bin_rhs = bin_edges[:-1], bin_edges[1:]
        bin_centers = (bin_lhs + bin_rhs) / 2

        shells.append(bin_centers[np.where(counts > 0)])

        nshells = [len(s) for s in shells]

        min_nshells = np.min(nshells)
        possible_shells = np.where(nshells == min_nshells)[0]
        chosen_shells = shells[possible_shells[-1]]

        # Force first shell to 0, round to nearest 100 otherwise
        chosen_shells[0] = 0
        chosen_shells = np.around(chosen_shells, -2).astype("int")

        # Get bval indices via mindist to shell
        rep_shells = np.tile(chosen_shells, [len(bvals), 1])
        rep_bvals = np.tile(bvals, [len(chosen_shells), 1]).T

        # Compute shell indices via abs diff between bvals and shells
        shell_ind = np.argmin(np.abs(rep_bvals - rep_shells), 1)
        shell_ind = chosen_shells[shell_ind]

        # Get list of indices shells to volumes
        shell_to_vol = dict()
        for shell in chosen_shells.tolist():
            shell_to_vol[shell] = np.where(shell_ind == int(shell))[0].tolist()

    return chosen_shells, shell_ind, shell_to_vol


def get_shells_from_bvals(bval_fname: str, out_json: str) -> None:
    bvals = np.loadtxt(bval_fname)
    out_dict = dict()
    # If single bval, set manually, otherwise find shells
    if bvals.size == 1:
        shells = [np.around(bvals, -2).astype("int").tolist()]
        out_dict["shells"] = shells
        out_dict["vol_to_shell"] = shells
        out_dict["shell_to_vol"] = {str(shells[0]): [0]}  # Point to index 0
    else:
        chosen_shells, shell_ind, shell_to_vol = find_shells(bvals)
        out_dict["shells"] = chosen_shells.tolist()
        out_dict["vol_to_shell"] = shell_ind.tolist()
        out_dict["shell_to_vol"] = shell_to_vol

    # Write shells and indices to json
    with open(out_json, "w") as out_fname:
        json.dump(out_dict, out_fname, indent=4)


if __name__ == "__main__":
    get_shells_from_bvals(
        bval_fname=snakemake.input.bval,  # noqa: F821
        out_json=snakemake.output.json,  # noqa: F821
    )
