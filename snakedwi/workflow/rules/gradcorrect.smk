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


def get_dwi_ref_for_gradcorrect(wildcards):
    if config["gradcorrect_coeffs"]:
        checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[
            0
        ]
        ([method],) = glob_wildcards(os.path.join(checkpoint_output, "sdc-{method}"))

        return sdc_methods[method]
    else:
        return ""


rule find_gradcorrect_warp:
    input:
        dwiref=get_dwi_ref_for_gradcorrect,
        grad_coeff=config["gradcorrect_coeffs"] or "",
    output:
        dwiref=bids(
            root=work,
            datatype="dwi",
            desc="gradcorrect",
            suffix="b0.nii.gz",
            **subj_wildcards,
        ),
        warp=bids(
            root=work,
            datatype="transforms",
            desc="gradcorrect",
            type_="fsl",
            suffix="xfm.nii.gz",
            **subj_wildcards,
        ),
    log:
        f"logs/find_gradcorrect_warp/{'.'.join(subj_wildcards.values())}.log",
    benchmark:
        f"benchmarks/find_gradcorrect_warp/{'.'.join(subj_wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    group:
        "subj"
    shadow:
        "minimal"
    params:
        fov=0.2,
        numpoints=120,
        out_warp="fullWarp_abs.nii.gz",
    container:
        config["singularity"]["gradcorrect"]
    shell:
        "gradient_unwarp.py {input.dwiref} {output.dwiref} siemens "
        "-g {input.grad_coeff} -n --fovmin -{params.fov} --fovmax {params.fov} "
        "--numpoints {params.numpoints} --verbose && "
        "mv {params.out_warp} {output.warp}"


rule get_jacobian_determinant:
    input:
        ref=rules.find_gradcorrect_warp.output.dwiref,
        warp=rules.find_gradcorrect_warp.output.warp,
    output:
        detjac=bids(
            root=root,
            datatype="dwi",
            desc="gradcorrect",
            suffix="detjac.nii.gz",
            **subj_wildcards,
        ),
    log:
        f"logs/get_jacobian_determinant/{'.'.join(subj_wildcards.values())}.log",
    benchmark:
        f"benchmarks/get_jacobian_determinant/{'.'.join(subj_wildcards.values())}.tsv"
    group:
        "subj"
    threads: 1
    shell:
        "reg_jacobian -ref {input.ref} -def {input.warp} -jac {output.detjac}"


rule convert_gradcorrect_to_itk:
    input:
        warp=rules.find_gradcorrect_warp.output.warp,
        refvol=rules.find_gradcorrect_warp.output.dwiref,
    output:
        warp=bids(
            root=root,
            datatype="transforms",
            desc="gradcorrect",
            type_="itk",
            suffix="xfm.nii.gz",
            **subj_wildcards,
        ),
    log:
        f"logs/convert_gradcorrect_to_itk/{'.'.join(subj_wildcards.values())}.log",
    benchmark:
        f"benchmarks/convert_gradcorrect_to_itk/{'.'.join(subj_wildcards.values())}.tsv"
    group:
        "subj"
    threads: 1
    container:
        config["singularity"]["connectome_workbench"]
    shell:
        "wb_command -convert-warpfield -from-fnirt {input.warp} {input.refvol} -absolute -to-itk {output.warp}"


rule gradcorrect_t1w:
    input:
        t1=bids(root=work, datatype="anat", **subj_wildcards, suffix="T1w.nii.gz"),
        grad_coeff=config["gradcorrect_coeffs"] or "",
    output:
        t1=bids(
            root=work,
            **subj_wildcards,
            desc="gradcorrect",
            datatype="anat",
            suffix="T1w.nii.gz",
        ),
        warp=bids(
            root=root,
            **subj_wildcards,
            suffix="xfm.nii.gz",
            desc="gradcorrect",
            space="T1w",
            datatype="transforms"
        ),
    log:
        f"logs/gradcorrect_t1w/{'.'.join(subj_wildcards.values())}.log",
    benchmark:
        f"benchmarks/gradcorrect_t1w/{'.'.join(subj_wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shadow:
        "minimal"
    params:
        fov=0.2,
        numpoints=120,
        out_warp="fullWarp_abs.nii.gz",
    container:
        config["singularity"]["gradcorrect"]
    group:
        "subj"
    shell:
        "gradient_unwarp.py {input.t1} {output.t1} siemens "
        "-g {input.grad_coeff} -n --fovmin -{params.fov} --fovmax {params.fov} "
        "--numpoints {params.numpoints} --verbose && "
        "mv {params.out_warp} {output.warp}"
