# this checkpoint allows us to use {shell} wildcard for downstream rules
checkpoint split_shell_avgs:
    input:
        avg_4d=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells.nii.gz",
            **subj_wildcards
        ),
    params:
        out_prefix=lambda wildcards, output: os.path.join(output.nii_dir, "dwi."),
    output:
        nii_dir=directory(
            bids(
                root="work",
                desc="topup",
                method="jac",
                suffix="dwi.avgshells",
                **subj_wildcards
            )
        ),
    container:
        config["singularity"]["fsl"]
    shell:
        "mkdir -p {output}; fslsplit {input} {params.out_prefix}"


# n4 on each shell
rule n4_shell_avg:
    input:
        avg_nii=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/dwi.{shell}.nii.gz",
            **subj_wildcards
        ),
    output:
        n4_nii=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/dwi_n4.{shell}.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["ants"]
    shell:
        "N4BiasFieldCorrection -i {input} -o {output}"


# rescale intensities, clip off first/last 5% of intensities, then rescale to 0-2000
rule rescale_shell_avg:
    input:
        in_nii="{infile}.{shell}.nii.gz",
    output:
        rescale_nii="{infile}_rescale.{shell}.nii.gz",
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d -verbose {input} -clip 5% 95% -stretch 0% 99% 0 2000 -o {output}"


# make initial brainmask much larger than needed (will be refined later.. )
rule bet_shell_avg:
    input:
        in_nii="{infile}.{shell}.nii.gz",
    params:
        frac=0.1,
    output:
        bet_nii="{infile}_bet.{shell}.nii.gz",
    container:
        config["singularity"]["fsl"]
    shell:
        "bet {input} {output} -f {params.frac}"


rule smooth_binarize_shell_avg:
    input:
        in_nii="{infile}.{shell}.nii.gz",
    params:
        smoothing="2mm",
        sdt_thresh=2,  # in mm, +ve for 2mm larger than bet mask
    output:
        mask_nii="{infile}_binarize.{shell}.nii.gz",
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input} -binarize -sdt -smooth {params.smoothing} -threshold {params.sdt_thresh} inf 0 1 -o {output}"


# now, run n4 using mask from b0 processing, with many more iterations
rule n4_shell_avg_withb0mask:
    input:
        avg_nii=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/dwi.{shell}.nii.gz",
            **subj_wildcards
        ),
        mask_nii=bids(
            root="work",
            suffix="mask.nii.gz",
            desc="brain",
            from_="avgb0",
            **subj_wildcards
        ),
    output:
        n4_nii=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/dwi_n4withb0mask.{shell}.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["ants"]
    shell:
        "N4BiasFieldCorrection -i {input.avg_nii} -o {output} -x {input.mask_nii} --convergence [200x200x200x200] "


def get_diffweighted_shells_for_tissue_seg(wildcards):
    checkpoint_output = checkpoints.split_shell_avgs.get(**wildcards).output[0]
    return sorted(
        expand(
            bids(
                root="work",
                desc="topup",
                method="jac",
                suffix="dwi.avgshells/dwi_n4_rescale.{shell}.nii.gz",
                **subj_wildcards
            ),
            **wildcards,
            shell=glob_wildcards(os.path.join(checkpoint_output, "dwi.{i}.nii.gz")).i[
                1:
            ]
        )
    )


# now, run tissue segmentation on all three shells using initial mask
rule tissue_seg_kmeans_init:
    input:
        shells=get_diffweighted_shells_for_tissue_seg,
        mask=bids(
            root="work",
            suffix="mask.nii.gz",
            desc="brain",
            from_="avgb0",
            **subj_wildcards
        ),
    params:
        intensity_images=lambda wildcards, input: " ".join(
            [f"--intensity-image [{img},1]" for img in input.shells]
        ),
        k="{k}",  # 4 should be enough to get bg as the first class
        posterior_fmt="posteriors_%d.nii.gz",
        posterior_glob="posteriors_*.nii.gz",
    output:
        seg=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_dseg.nii.gz",
            **subj_wildcards
        ),
        posteriors=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_probseg.nii.gz",
            **subj_wildcards
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["prepdwi"]
    shell:
        #merge posteriors into a 4d file (intermediate files will be removed b/c shadow) - TODO: update this so just requires ants
        "Atropos -d 3  {params.intensity_images} -i KMeans[{params.k}] -x {input.mask} -o [{output.seg},{params.posterior_fmt}] && "
        "fslmerge -t {output.posteriors} {params.posterior_glob} "


# get class 0 from tissue seg - corresponds to lowest intensity
rule extract_posterior_bgnd:
    input:
        posterior_4d=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_probseg.nii.gz",
            **subj_wildcards
        ),
    output:
        posterior_bgnd=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_label-bg_probseg.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["fsl"]
    shell:
        "fslroi {input} {output} 0 1"


# subtract posterior probability of tissue class 1 from the initial brainmask
rule refine_mask_with_tissue_prob:
    input:
        mask=bids(
            root="work",
            suffix="mask.nii.gz",
            desc="brain",
            from_="avgb0",
            **subj_wildcards
        ),
        posterior_bgnd=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_label-bg_probseg.nii.gz",
            **subj_wildcards
        ),
    output:
        mask=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_label-brain_probseg.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["fsl"]
    shell:
        "fslmaths {input.mask} -sub {input.posterior_bgnd} {output.mask}"


# smooth and threshold to binarize
rule smooth_threshold_refined_mask:
    input:
        mask=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_label-brain_probseg.nii.gz",
            **subj_wildcards
        ),
    params:
        smooth="{smooth}",
    output:
        mask=bids(
            root="work",
            desc="topup",
            method="jac",
            suffix="dwi.avgshells/atropos_k-{k}_initmasking_label-brain_smooth-{smooth}_mask.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.mask} -smooth {params.smooth} -threshold 0.5 Inf 1 0 -o {output.mask}"
