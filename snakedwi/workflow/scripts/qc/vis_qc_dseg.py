import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
import skimage as skimg
from nilearn import plotting

matplotlib.use("Agg")


with plt.ioff():
    mask = nib.load(snakemake.input["seg"])
    background = nib.load(snakemake.input["img"])
    boundary = nib.Nifti1Image(
        skimg.segmentation.find_boundaries(
            mask.get_fdata(),
            mode="inner",
        ).astype("float32"),
        mask.affine,
    )

    # Static Edge overlay
    ax_size = 2
    n_slices = 7
    fig = plt.figure(
        figsize=(ax_size * n_slices, ax_size * 4),
        layout="constrained",
        facecolor="k",
        edgecolor="k",
    )
    ax = fig.add_subplot(111)
    display = plotting.plot_img(
        background,
        n_slices,
        display_mode="mosaic",
        black_bg=True,
        colorbar=False,
        annotate=False,
        cmap="gray",
        axes=ax,
    )
    display.add_overlay(boundary, cmap="autumn")
    ax.axis("off")
    fig.savefig(snakemake.output.png, dpi=600)
    plt.close(fig)

    html_view = plotting.view_img(
        stat_map_img=boundary,
        bg_img=background,
        opacity=0.5,
        cmap="autumn",
        threshold=0.5,
        colorbar=False,
        symmetric_cmap=False,
        resampling_interpolation="nearest",
        title="sub-{subject}".format(**snakemake.wildcards),
    )

    html_view.save_as_html(snakemake.output.html)
