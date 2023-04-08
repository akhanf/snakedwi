#!/usr/bin/env python

from typing import List

import matplotlib
from nilearn import plotting


def vis_regqc(
    anat_nii: str,
    ref: str,
    out_html: str,
    out_png: str,
    smk_wildcards: List(str),
) -> None:

    matplotlib.use("Agg")

    # HTML Output
    html_view = plotting.view_img(
        stat_map_img=ref,
        bg_img=anat_nii,
        opacity=0.5,
        cmap="viridis",
        dim=-1,
        symmetric_cmap=False,
        title="sub-{subject}".format(**smk_wildcards),
    )
    html_view.save_as_html(out_html)

    # PNG Output
    display = plotting.plot_anat(anat_nii, display_mode="ortho")
    display.add_contours(ref, colors="r")
    display.savefig(out_png)
    display.close()


if __name__ == "__main__":
    vis_regqc(
        anat_nii=snakemake.input.flo,  # noqa: F821
        ref=snakemake.input.ref,  # noqa: F821
        out_html=snakemake.output.html,  # noqa: F821
        out_png=snakemake.output.png,  # noqa: F821
        smk_wildcards=snakemake.wildcards,  # noqa: F821
    )
