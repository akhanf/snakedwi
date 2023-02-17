
rule find_gradcorrect_warp:
    input:
        data=rules.cp_dwi_ref.output,
        grad_coeff=config["gradcorrect_coeffs"] or "",
    output:
        data=bids(
            root=work,
            datatype="dwi",
            desc="unwarped",
            suffix="b0.nii.gz",
            **subj_wildcards,
        ),
        warp=bids(
            root=work,
            datatype="dwi",
            desc="gradCorrect",
            type_="fsl",
            suffix="warp.nii.gz",
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
        "gradient_unwarp.py {input.data} {output.data} siemens "
        "-g {input.grad_coeff} -n --fovmin -{params.fov} --fovmax {params.fov} "
        "--numpoints {params.numpoints} --verbose && "
        "mv {params.out_warp} {output.warp}"



rule get_jacobian_determinant:
    input:
        ref=rules.find_gradcorrect_warp.output.data,
        warp=rules.find_gradcorrect_warp.output.warp,
    output:
        detjac=bids(
            root=root,
            datatype="dwi",
            desc="gradCorrect",
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
        refvol=rules.cp_dwi_ref.output,
    output:
        warp=bids(
            root=work,
            datatype="dwi",
            desc="gradCorrect",
            type_="itk",
            suffix="warp.nii.gz",
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


rule resample_dwi_to_t1w:
    input:
        ref=bids(
            root=work,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="dwiref",
            proc="crop",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="eddy",
            datatype="dwi",
            **subj_wildcards
        ),
        xfm_itk=bids(
            root=work,
            suffix="xfm.txt",
            from_="dwi",
            to="T1w",
            type_="itk",
            datatype="dwi",
            **subj_wildcards
        ),
        **(
            {"gradcorrect_warp": rules.convert_gradcorrect_to_itk.output}
            if config["gradcorrect_coeffs"] else {}
        ),
    params:
        interpolation="Linear",
    output:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="eddy",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["ants"]
    resources:
        mem_mb=32000,  #-- this is going to be dependent on size of image.. 
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 --input-image-type 3 "
        "--input {input.dwi} --reference-image {input.ref} "
        "--transform {input.xfm_itk} "
        f"{'-t {input.gradcorrect_warp}' if config['gradcorrect_coeffs'] else ''} "
        "--interpolation {params.interpolation} --output {output.dwi} --verbose "

rule gradcorrect_t1w:
    input:
        t1=rules.import_t1.output,
        grad_coeff=config["gradcorrect_coeffs"] or "",
    output:
        t1=bids(
            root=root,
            **subj_wildcards,
            desc="preproc",
            datatype="anat",
            suffix="T1w.nii.gz",
        ),
        warp=bids(
            root=root,
            **subj_wildcards,
            suffix="warp.nii.gz",
            desc="gradCorrect",
            space="T1w",
            datatype="anat"
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


