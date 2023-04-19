from nilearn import plotting
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib


matplotlib.use("Agg")


with plt.ioff():
    lines = nib.load(snakemake.input["lines"])
    overlay = nib.load(snakemake.input["overlay"])
    background = nib.load(snakemake.input["background"])

    # make sure sform is set (seems to cause problems with nilearn if unset)
    for img in (lines, overlay, background):
        img.set_sform(img.affine)


    # Static Edge overlay
    ax_size = 2
    n_slices = 7
    fig = plt.figure(
        figsize=(ax_size*n_slices, ax_size*4),
        layout="constrained",
        facecolor='k',
        edgecolor='k'
    )
    ax = fig.add_subplot(111)
    display = plotting.plot_img(
        overlay,
        n_slices,
        display_mode='mosaic',
        black_bg=True,
        colorbar=False,
        annotate=False,
        cmap='gray',
        axes=ax,
    )
    display.add_overlay(lines, cmap='autumn')
    ax.axis('off')
    fig.savefig(snakemake.output.png, dpi=600)
    plt.close(fig)


    # Interactive heatmap overlay
    html_view = plotting.view_img(
        stat_map_img=overlay,
        bg_img=background,
        opacity=0.5,
        cmap="hot",
        threshold="1%",
        symmetric_cmap=False,
        title="sub-{subject}".format(**snakemake.wildcards),
    )
    html_view.save_as_html(snakemake.output.html)
        