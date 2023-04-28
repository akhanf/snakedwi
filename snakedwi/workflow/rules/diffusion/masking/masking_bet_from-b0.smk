# n4
rule n4_avg_b0:
    input:
        dwiref=rules.cp_dwi_ref.output.dwi_ref,
    output:
        n4_avgb0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="n4",
            datatype="dwi",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "N4BiasFieldCorrection -i {input} -o {output}"


# rescale intensities - clip first/last 5% of intensities -> rescale to 0-2000
rule rescale_avg_b0:
    input:
        n4_avgb0=rules.n4_avg_b0.output.n4_avgb0,
    output:
        rescale_b0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="rescale",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d -verbose {input} -clip 5% 95% -stretch 0% 99% 0 2000 -o {output}"


rule bet_avg_b0:
    input:
        rescale_b0=rules.rescale_avg_b0.output.rescale_b0,
    params:
        bet_frac=config["b0_bet_frac"],
    output:
        b0_brain=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="bet",
            datatype="dwi",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["fsl"]
    group:
        "subj"
    shell:
        "bet {input} {output} -f {params.bet_frac}"


rule binarize_avg_b0:
    input:
        rules.bet_avg_b0.output.b0_brain,
    output:
        mask=bids(
            root=work,
            suffix="mask.nii.gz",
            desc="brain",
            method="bet_from-b0",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -binarize  -o {output}"
