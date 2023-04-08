rule run_synthSR:
    input:
        "{prefix}.nii.gz",
    output:
        "{prefix}SynthSR.nii.gz",
    threads: 8
    group:
        "subj"
    container:
        config["singularity"]["synthsr"]
    shadow:
        "minimal"
    shell:
        "python /SynthSR/scripts/predict_command_line.py "
        "--cpu --threads {threads} "
        "{input} {output}"


rule reslice_synthSR_b0:
    input:
        ref=rules.moco_subj_bzeros_4d.output.nii_avg3d,
        synthsr=bids(
            root=work,
            suffix="b0SynthSR.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
    output:
        synthsr=bids(
            root=work,
            suffix="b0SynthSRresliced.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input.ref} {input.synthsr} -reslice-identity -o {output.synthsr}"


rule rigid_reg_t1_to_b0_synthsr:
    input:
        flo=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1wSynthSR.nii.gz" ** subj_wildcards,
        ),
        ref=rules.reslice_synthSR_b0.output.synthsr,
    output:
        warped_subj=bids(
            root=work,
            suffix="T1wSynthSRreg.nii.gz",
            space="rigidb0",
            datatype="dwi",
            **subj_wildcards
        ),
        xfm_ras=bids(
            root=work,
            suffix="xfm.txt",
            from_="T1wSynthSR",
            to="b0SynthSR",
            type_="ras",
            desc="rigid",
            datatype="transforms",
            **subj_wildcards
        ),
    container:
        config["singularity"]["prepdwi"]
    log:
        bids(
            root="logs",
            suffix="rigid_reg_t1_to_b0_synthsr.txt",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} "
        "-res {output.warped_subj} -aff {output.xfm_ras} > {log}"


def get_restrict_deformation(wildcards, input, output):
    peaxis = get_pe_axis(wildcards)
    if peaxis == "i":
        return "1x0x0"
    elif peaxis == "j":
        return "0x1x0"
    elif peaxis == "k":
        return "0x0x1"  # unlikely..


rule reg_b0_to_t1_synthsr:
    input:
        rules.check_subj_dwi_metadata.output,
        t1synth=rules.rigid_reg_t1_to_b0_synthsr.output.warped_subj,
        b0synth=rules.reslice_synthSR_b0.output.synthsr,
    params:
        general_opts="-d 3 -v",
        metric=lambda wildcards, input: (
            f"MeanSquares[{input.t1synth},{input.b0synth}]"
        ),
        transform="SyN[0.1]",
        convergence="100x50x20",
        smoothing_sigmas="2x1x0vox",
        shrink_factors="4x2x1",
        restrict_deformation=get_restrict_deformation,
    # "0x1x0", should be set to the phase encode dir - can read json for this.
    output:
        unwarped=bids(
            root=work,
            suffix="b0SynthSRunwarped.nii.gz",
            desc="Syn",
            datatype="dwi",
            **subj_wildcards
        ),
        fwd_xfm=bids(
            root=root,
            suffix="xfm.nii.gz",
            from_="b0SynthSR",
            to="T1wSynthSR",
            type_="itk",
            desc="SyN",
            datatype="transforms",
            **subj_wildcards
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["ants"]
    log:
        bids(root="logs", suffix="reg_b0_to_t1_synthsr.txt", **subj_wildcards),
    group:
        "subj"
    shell:
        "antsRegistration {params.general_opts} --metric {params.metric} "
        "--convergence {params.convergence} "
        "--smoothing-sigmas {params.smoothing_sigmas} "
        "--shrink-factors {params.shrink_factors} "
        "--transform {params.transform} "
        "--output [ants_,{output.unwarped}] "
        "--restrict-deformation {params.restrict_deformation} > {log} && "
        "mv ants_0Warp.nii.gz {output.fwd_xfm} "


rule displacement_field_to_fmap:
    input:
        disp_field=rules.reg_b0_to_t1_synthsr.output.fwd_xfm,
        dwi_nii=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.nii.gz",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        )[0],
        dwi_json=lambda wildcards: get_dwi_indices(
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
    params:
        demean=True,
    output:
        fmap=bids(
            root=work,
            datatype="dwi",
            suffix="fmap.nii.gz",
            desc="b0",
            method="synthSRsdc",
            **subj_wildcards
        ),
    container:
        config["singularity"]["python"]
    script:
        "../scripts/diffusion/sdc/displacement_field_to_fmap.py"


rule apply_unwarp_synthsr:
    input:
        b0=rules.moco_subj_bzeros_4d.output.nii_avg3d,
        xfm=rules.reg_b0_to_t1_synthsr.output.fwd_xfm,
    output:
        unwarped=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="unwarped",
            method="synthSRsdc",
            **subj_wildcards
        ),
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 -i {input.b0} -o {output.unwarped} "
        "-r {input.b0} -t {input.xfm} -e 0 -n NearestNeighbor"
