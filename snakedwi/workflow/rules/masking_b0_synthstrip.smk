rule synthstrip_b0:
    input:
        b0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        temp(
            bids(
                root=work,
                suffix="mask.nii.gz",
                desc="nofixhdrbrain",
                method="synthstrip",
                from_="{desc}b0",
                datatype="dwi",
                **subj_wildcards
            )
        ),
    container:
        config["singularity"]["synthstrip"]
    threads: 8
    group:
        "subj"
    shell:
        "python3 /freesurfer/mri_synthstrip -i {input} -m {output}"


rule synthstrip_b0_fix_header:
    input:
        b0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
        mask=bids(
            root=work,
            suffix="mask.nii.gz",
            desc="nofixhdrbrain",
            method="synthstrip",
            from_="{desc}b0",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        mask=bids(
            root=work,
            suffix="mask.nii.gz",
            desc="brain",
            method="synthstrip",
            from_="{desc}b0",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.b0} {input.mask} -copy-transform -o {output.mask}"
