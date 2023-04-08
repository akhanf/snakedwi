# just grab the first T1w for now:
rule import_t1:
    input:
        nii=lambda wildcards: expand(
            input_path["T1w"],
            zip,
            **filter_list(input_zip_lists["T1w"], wildcards)
        )[0],
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input.nii} {output.nii}"


def get_input_for_synthstrip(wildcards):
    if config["gradcorrect_coeffs"]:
        return bids(
            root=work,
            datatype="anat",
            desc="gradcorrect",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        )
    else:
        return bids(
            root=work,
            datatype="anat",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        )


rule synthstrip_t1:
    input:
        t1=get_input_for_synthstrip,
    output:
        mask=temp(
            bids(
                root=work,
                datatype="anat",
                desc="nofixhdrbrain",
                suffix="mask.nii.gz",
                **subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["synthstrip"]
    threads: 8
    shadow:
        "minimal"
    shell:
        "python3 /freesurfer/mri_synthstrip -i {input.t1} -m {output.mask}"


rule fixheader_synthstrip:
    input:
        t1=rules.import_t1.output.nii,
        mask=rules.synthstrip_t1.output.mask,
    output:
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.t1} {input.mask} -copy-transform -o {output.mask}"


rule n4_t1_withmask:
    input:
        t1=rules.import_t1.output.nii,
        mask=rules.fixheader_synthstrip.output.mask,
    output:
        t1=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
    threads: 8
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}"
