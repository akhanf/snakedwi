sdc_methods = {
    "topup": bids(
        root=work,
        suffix="b0.nii.gz",
        desc="topup",
        method="jac",
        datatype="dwi",
        **subj_wildcards,
    ),
    "synthsr": bids(
        root=work,
        datatype="dwi",
        suffix="b0.nii.gz",
        desc="unwarped",
        method="synthSRsdc",
        **subj_wildcards,
    ),
    "sdcflow": bids(
        root=work,
        datatype="dwi",
        suffix="b0.nii.gz",
        desc="unwarped",
        method="synsdc",
        **subj_wildcards,
    ),
    "none": bids(
        root=work,
        suffix="b0.nii.gz",
        datatype="dwi",
        desc="moco",
        **subj_wildcards,
    ),
}


def get_dwi_ref(wildcards):
    method = get_sdc_method(wildcards)

    if config["gradcorrect_coeffs"]:
        # if using gradcorrect, checkpoint for sdc method will be
        # done in the find_gradcorrect_warp rule instead
        return bids(
            root=work,
            datatype="dwi",
            desc="gradcorrect",
            suffix="b0.nii.gz",
            **subj_wildcards,
        )
    else:
        return sdc_methods[method]


def get_dwi_ref(wildcards):
    if config["gradcorrect_coeffs"]:
        # if using gradcorrect, then the checkpoint for sdc method will be
        # done in the find_gradcorrect_warp rule instead
        return bids(
            root=work,
            datatype="dwi",
            desc="gradcorrect",
            suffix="b0.nii.gz",
            **subj_wildcards,
        )
    else:
        checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[
            0
        ]
        ([method],) = glob_wildcards(os.path.join(checkpoint_output, "sdc-{method}"))

        return sdc_methods[method]


rule cp_dwi_ref:
    input:
        get_dwi_ref,
    output:
        dwi_ref=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="dwiref",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"
