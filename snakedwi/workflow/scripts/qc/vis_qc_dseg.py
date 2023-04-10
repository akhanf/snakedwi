#!/usr/bin/env python
import matplotlib
from nilearn import plotting
from numpy import ndarray


def vis_qc_dseg(
    dseg: str, bg_img: ndarray, out_html: str, out_png: str, wildcards: dict
) -> None:
    matplotlib.use("Agg")

    # Create HTML
    html_view = plotting.view_img(
        stat_map_img=dseg,
        bg_img=bg_img,
        opacity=0.5,
        cmap="viridis",
        dim=-1,
        threshold=0.5,
        symmetric_cmap=False,
        title="sub-{subject}".format(**wildcards),
    )
    html_view.save_as_html(out_html)

    # Create PNG
    display = plotting.plot_roi(
        roi_img=dseg, bg_img=bg_img, display_mode="ortho"
    )
    display.savefig(out_png)
    display.close()


if __name__ == "__main__":
    vis_qc_dseg(
        dseg=snakemake.input.seg,  # noqa: F821
        bg_img=snakemake.input.img,  # noqa: F821
        out_html=snakemake.output.html,  # noqa: F821
        out_png=snakemake.output.png,  # noqa: F821
        wildcards=snakemake.wildcards,  # noqa: F821
    )
