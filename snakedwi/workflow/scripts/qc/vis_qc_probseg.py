#!/usr/bin/env python

import matplotlib
from nilearn import plotting


def vis_qc_probseg(img: str, seg4d: str, out_png: str) -> None:
    matplotlib.use("Agg")

    # Save output png
    display = plotting.plot_prob_atlas(bg_img=img, maps_img=seg4d, display_mode="ortho")
    display.savefig(out_png)
    display.close()


if __name__ == "__main__":
    vis_qc_probseg(
        img=snakemake.input.img,  # noqa: F821
        seg4d=snakemake.input.seg4d,  # noqa: F821
        out_png=snakemake.output.png,  # noqa: F821
    )
