rule extract_phenc:
    input:
        phenc_concat=rules.concat_phase_encode_txt.output.phenc_concat
    output:
        synb0_phenc=bids(
            root=work,
            datatype="dwi",
            desc="synb0",
            suffix="phenc.txt",
            **subj_wildcards
        )
    container:
        config["singularity"]["python"]
    group:
        "subj"
    script:
        "../../../scripts/diffusion/sdc/get_synb0_phase_encode.py"


# Symlinks don't work for synb0
rule setup_synb0:
    input:
        synb0_phenc=rules.extract_phenc.output.synb0_phenc,
        b0=rules.moco_subj_bzeros_4d.output.nii_avg3d,
        t1=rules.n4_t1_withmask.output.t1
    output:
        synb0_indir=directory(
            Path(
                bids(
                    root=work,
                    datatype="synb0_in",
                    **subj_wildcards
                )
            ).parent
        )
    group:
        "subj"
    shell:
        "mkdir -p {output.synb0_indir} && "
        "cp {input.b0} {output.synb0_indir}/b0.nii.gz && "
        "cp {input.t1} {output.synb0_indir}/T1.nii.gz && "
        "cp {input.synb0_phenc} {output.synb0_indir}/acqparams.txt"


rule run_synb0:
    input:
        synb0_indir=rules.setup_synb0.output.synb0_indir,
    params:
        fs_license=config["fs_license"],
        container=config["singularity"]["synb0"],
        out_dir=directory(
            Path(
                bids(
                    root=work,
                    datatype="synb0_out",
                    **subj_wildcards,
                )
            ).parent
        )
    output:
        b0_all=bids(
            root=work,
            suffix="b0all.nii.gz",
            desc="synb0topup",
            datatype="dwi",
            **subj_wildcards
        ),
        fieldcoef=bids(
            root=work,
            desc="topup",
            suffix="fieldcoef.nii.gz",
            method="synb0",
            datatype="dwi",
            **subj_wildcards
        ),
        movpar=bids(
            root=work,
            desc="topup",
            suffix="movpar.txt",
            method="synb0",
            datatype="dwi",
            **subj_wildcards
        ),
    log:
        bids(root="logs", suffix="synb0.log", **subj_wildcards),
    group: 
        "subj"
    shell:
        "mkdir -p {params.out_dir} && "
        "singularity run -e -B {input.synb0_indir}/:/INPUTS "
        "-B {params.out_dir}/:/OUTPUTS "
        "-B {params.fs_license}:/extra/freesurfer/license.txt "
        "{params.container} --stripped 2> {log} && "
        "mv {params.out_dir}/b0_all_topup.nii.gz {output.b0_all} && "
        "mv {params.out_dir}/topup_fieldcoef.nii.gz {output.fieldcoef} && "
        "mv {params.out_dir}/topup_movpar.txt {output.movpar}"


# Extract unwarped subject (non-syn) b0
rule get_ref_synb0:
    input:
        b0_all=rules.run_synb0.output.b0_all
    output:
        b0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="topup",
            method="synb0",
            datatype="dwi",
            **subj_wildcards,
        )
    container:
        config["singularity"]["mrtrix"]
    group:
        "subj"
    shell:
        "mrconvert {input.b0_all} -coord 3 0 -axes 0,1,2 {output.b0}"