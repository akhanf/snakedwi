rule reg_dwi_to_t1:
    input:
        t1wsynth=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1wSynthSR.nii.gz",
            **subj_wildcards
        ),
        avgb0synth=bids(
            root=work,
            datatype="dwi",
            desc="dwiref",
            suffix="b0SynthSR.nii.gz",
            **subj_wildcards
        ),
        avgb0=rules.cp_dwi_ref.output.dwi_ref,
    params:
        general_opts="-d 3",
        rigid_opts=(
            "-m NMI -a -dof 6 -ia-{rigid_dwi_t1_init} -n {rigid_dwi_t1_iters}"
        ).format(
            rigid_dwi_t1_init=config["rigid_dwi_t1_init"],
            rigid_dwi_t1_iters=config["rigid_dwi_t1_iters"],
        ),
    output:
        warped_avgb0=bids(
            root=work,
            datatype="dwi",
            space="T1w",
            desc="dwiref",
            suffix="avgb0.nii.gz",
            **subj_wildcards
        ),
        xfm_ras=bids(
            root=work,
            datatype="transforms",
            from_="dwi",
            to="T1w",
            type_="ras",
            suffix="xfm.txt",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    log:
        bids(
            root="logs",
            suffix="reg_b0_to_t1.txt",
            datatype="dwi",
            **subj_wildcards
        ),
    threads: 8
    shell:
        """
        greedy -threads {threads} {params.general_opts} {params.rigid_opts} \\
            -i {input.t1wsynth} {input.avgb0synth} -o {output.xfm_ras} \\
            &> {log}

        greedy -threads {threads} {params.general_opts} -rf {input.t1wsynth} \\
            -rm {input.avgb0} {output.warped_avgb0} -r {output.xfm_ras} \\
             &>> {log}
        """


# rule qc_reg_dwi_t1:
#     input:
#         ref=bids(
#             root=work,
#             suffix="avgb0.nii.gz",
#             space="T1w",
#             desc="dwiref",
#             datatype="dwi",
#             **subj_wildcards
#         ),
#         flo=bids(
#             root=root,
#             suffix="T1w.nii.gz",
#             desc="preproc",
#             datatype="anat",
#             **subj_wildcards
#         ),
#     output:
#         png=report(
#             bids(
#                 root=root,
#                 datatype="qc",
#                 suffix="reg.png",
#                 **subj_wildcards,
#                 from_="dwiref",
#                 to="T1w"
#             ),
#             caption="../report/reg_dwi_t1.rst",
#             category="B0 T1w registration",
#         ),
#         html=bids(
#             root=root,
#             datatype="qc",
#             suffix="reg.html",
#             from_="dwiref",
#             to="T1w",
#             **subj_wildcards
#         ),
#     group:
#         "subj"
#     container:
#         config["singularity"]["python"]
#     script:
#         "../scripts/vis_regqc.py"

rule qc_get_t1_edges:
    input:
        brain=rules.n4_t1_withmask.output.t1,
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            **subj_wildcards,
        )
    output:
        temp(
            bids(
                work,
                datatype="anat",
                desc="cannyEdge",
                suffix="T1w.nii.gz",
                **subj_wildcards,
            )
        )
    group: 'subj'
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input} -times -canny 1mm 10 20 {output}"


rule qc_reg_dwi_t1:
    input:
        overlay=rules.reg_dwi_to_t1.output.warped_avgb0,
        background=rules.n4_t1_withmask.output.t1,
        lines=rules.qc_get_t1_edges.output[0],
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                from_="dwiref",
                to="T1w",
                suffix="reg.png",
                **subj_wildcards,
            ),
            caption="../report/reg_dwi_t1.rst",
            category="B0 T1w registration",
        ),
        html=bids(
            root=root,
            datatype="qc",
            from_="dwiref",
            to="T1w",
            suffix="reg.html",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/qc/vis_regqc.py"



rule convert_xfm_ras2itk:
    input:
        xfm_ras=rules.reg_dwi_to_t1.output.xfm_ras,
    output:
        xfm_itk=bids(
            root=root,
            datatype="transforms",
            from_="dwi",
            to="T1w",
            type_="itk",
            suffix="xfm.txt",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule convert_xfm_ras2fsl:
    input:
        t1w=rules.n4_t1_withmask.output.t1,
        avgb0=rules.cp_dwi_ref.output.dwi_ref,
        xfm_ras=rules.reg_dwi_to_t1.output.xfm_ras,
    output:
        xfm_fsl=bids(
            root=work,
            datatype="transforms",
            from_="dwi",
            to="T1w",
            type_="fsl",
            suffix="xfm.txt",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input.xfm_ras} -ref {input.t1w} -src {input.avgb0} "
        "-ras2fsl -o {output.xfm_fsl}"


# Tight crop around b0 after rotating into T1w space
rule create_cropped_ref:
    input:
        warped_avgb0=rules.reg_dwi_to_t1.output.warped_avgb0,
    output:
        cropped_avgb0=bids(
            root=work,
            datatype="dwi",
            space="T1w",
            desc="dwiref",
            proc="crop",
            suffix="avgb0.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -trim 0vox {output}"


# for later resampling
rule write_nii_resolution_to_txt:
    input:
        nii="{prefix}.nii.gz",
    output:
        txt_fname="{prefix}.resolution_mm.txt",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/metadata/write_nii_resolution_to_txt.py"


# rules for creating reference image for each resampling scheme
# (only the rules that are required will be run)
rule create_cropped_ref_t1_resolution:
    input:
        cropped_avgb0=rules.create_cropped_ref.output.cropped_avgb0,
    output:
        avgb0_crop_resample=bids(
            root=work,
            datatype="dwi",
            space="T1w",
            desc="dwiref",
            proc="crop",
            res="T1w",
            suffix="avgb0.nii.gz",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule create_cropped_ref_dwi_resolution:
    input:
        cropped=(
            rules.create_cropped_ref_t1_resolution.output.avgb0_crop_resample
        ),
        res_txt=bids(
            root=work,
            datatype="dwi",
            desc="dwiref",
            suffix="b0.resolution_mm.txt",
            **subj_wildcards
        ),
    output:
        resampled=bids(
            root=work,
            datatype="dwi",
            space="T1w",
            desc="dwiref",
            proc="crop",
            res="orig",
            suffix="avgb0.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input.cropped} -resample-mm `cat {input.res_txt}` {output}"


rule create_cropped_ref_custom_resolution:
    input:
        cropped=(
            rules.create_cropped_ref_t1_resolution.output.avgb0_crop_resample
        ),
    params:
        resolution="x".join(
            [str(vox) for vox in config["resample_dwi"]["resample_mm"]]
        )
        + "mm",
    output:
        resampled=bids(
            root=work,
            datatype="dwi",
            space="T1w",
            desc="dwiref",
            proc="crop",
            res="custom",
            suffix="avgb0.nii.gz",
            **subj_wildcards
        ),
    container:
        config["singularity"]["itksnap"]
    group:
        "subj"
    shell:
        "c3d {input} -resample-mm {params.resolution} {output}"


rule resample_dwi_to_t1w:
    input:
        **(
            {"gradcorrect_warp": rules.convert_gradcorrect_to_itk.output}
            if config["gradcorrect_coeffs"]
            else {}
        ),
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
        xfm_itk=rules.convert_xfm_ras2itk.output.xfm_itk,
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
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    shell:
        "antsApplyTransforms -d 3 --input-image-type 3 "
        "--input {input.dwi} --reference-image {input.ref} "
        "--transform {input.xfm_itk} "
        f"{'-t {input.gradcorrect_warp}' if config['gradcorrect_coeffs'] else ''} "
        "--interpolation {params.interpolation} --output {output.dwi} --verbose "


rule resample_brainmask_to_t1w:
    input:
        **(
            {"gradcorrect_warp": rules.convert_gradcorrect_to_itk.output}
            if config["gradcorrect_coeffs"]
            else {}
        ),
        ref=rules.resample_dwi_to_t1w.input.ref,
        brainmask=get_b0_mask(),
        xfm_itk=rules.convert_xfm_ras2itk.output.xfm_itk,
    params:
        interpolation="NearestNeighbor",
    output:
        brainmask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
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
        "antsApplyTransforms -d 3 --input-image-type 0 "
        "--input {input.brainmask} --reference-image {input.ref} "
        "--transform {input.xfm_itk} "
        f"{'-t {input.gradcorrect_warp}' if config['gradcorrect_coeffs'] else ''} "
        "--interpolation {params.interpolation} --output {output.brainmask} --verbose "


rule rotate_bvecs_to_t1w:
    input:
        bvecs=bids(
            root=root,
            suffix="dwi.bvec",
            desc="eddy",
            datatype="dwi",
            **subj_wildcards
        ),
        bvals=bids(
            root=root,
            suffix="dwi.bval",
            desc="eddy",
            datatype="dwi",
            **subj_wildcards
        ),
        xfm_fsl=rules.convert_xfm_ras2fsl.output.xfm_fsl,
        script=os.path.join(
            workflow.basedir, f"scripts/diffusion/rotate_bvecs.sh"
        ),
    output:
        bvecs=bids(
            root=root,
            suffix="dwi.bvec",
            desc="eddy",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bvals=bids(
            root=root,
            suffix="dwi.bval",
            desc="eddy",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["fsl"]
    group:
        "subj"
    shell:
        "chmod a+x {input.script} && "
        "{input.script} {input.bvecs} {input.xfm_fsl} {output.bvecs} && "
        "cp -v {input.bvals} {output.bvals}"


# dti fitting on dwi in t1w space
rule dtifit_resampled_t1w:
    input:
        dwi=rules.resample_dwi_to_t1w.output.dwi,
        bvals=rules.rotate_bvecs_to_t1w.output.bvals,
        bvecs=rules.rotate_bvecs_to_t1w.output.bvecs,
        brainmask=rules.resample_brainmask_to_t1w.output.brainmask,
    params:
        out_basename=lambda wildcards, output: os.path.join(
            output.out_folder, "dti"
        ),
    output:
        out_folder=directory(
            bids(
                root=root,
                suffix="dtifit",
                desc="eddy",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            )
        ),
        out_fa=os.path.join(
            directory(
                bids(
                    root=root,
                    suffix="dtifit",
                    desc="eddy",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
            "dti_FA.nii.gz",
        ),
    container:
        config["singularity"]["fsl"]
    group:
        "subj"
    shell:
        "mkdir -p {output.out_folder} && "
        "dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} "
        "--mask={input.brainmask} --out={params.out_basename}"
