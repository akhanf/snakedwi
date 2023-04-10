def get_b0_init_mask():
    return bids(
        root=work,
        suffix="mask.nii.gz",
        desc="brain",
        method="synthstrip",
        from_="mocob0",
        datatype="dwi",
        **subj_wildcards,
    )


rule syn_sdc:
    input:
        in_epis=rules.moco_subj_bzeros_4d.output.nii_avg3d,
        in_anat=rules.n4_t1_withmask.output.t1,
        mask_anat=rules.fixheader_synthstrip.output.mask,
        std2anat_xfm=rules.invert_subj_to_template_xfm.output.xfm,
        epi_mask=get_b0_init_mask(),
        json_file=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.json",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        )[0],
    output:
        unwarped=bids(
            root=work,
            datatype="dwi",
            suffix="b0.nii.gz",
            desc="unwarped",
            method="synsdc",
            **subj_wildcards
        ),
        unwarped_mask=bids(
            root=work,
            datatype="dwi",
            suffix="mask.nii.gz",
            desc="brain",
            method="synsdc",
            **subj_wildcards
        ),
        xfm=bids(
            root=work,
            datatype="dwi",
            suffix="xfm.nii.gz",
            desc="itk",
            method="synsdc",
            **subj_wildcards
        ),
        fmap=bids(
            root=work,
            datatype="dwi",
            suffix="fmap.nii.gz",
            desc="b0",
            method="synsdc",
            **subj_wildcards
        ),
    threads: 8
    group:
        "subj"
    shadow:
        "minimal"
    container:
        config["singularity"]["sdcflows"]
    script:
        "../../../scripts/diffusion/sdc/sdcflows_syn.py"
