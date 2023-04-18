rule copy_inputs_for_bedpost:
    input:
        dwi=rules.resample_dwi_to_t1w.output.dwi,
        bval=rules.rotate_bvecs_to_t1w.output.bvals,
        bvec=rules.rotate_bvecs_to_t1w.output.bvecs,
        brainmask=rules.resample_brainmask_to_t1w.output.brainmask,
    output:
        diff_dir=temp(
            directory(
                bids(
                    root=work,
                    desc="eddy",
                    suffix="diffusion",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            )
        ),
        dwi=os.path.join(
            bids(
                root=work,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            "data.nii.gz",
        ),
        brainmask=os.path.join(
            bids(
                root=work,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            "nodif_brain_mask.nii.gz",
        ),
        bval=os.path.join(
            bids(
                root=work,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            "bvals",
        ),
        bvec=os.path.join(
            bids(
                root=work,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            "bvecs",
        ),
    group:
        "subj"
    shell:
        "mkdir -p {output.diff_dir} && "
        "cp {input.dwi} {output.dwi} && "
        "cp {input.brainmask} {output.brainmask} && "
        "cp {input.bval} {output.bval} && "
        "cp {input.bvec} {output.bvec} "


if config["use_bedpost_gpu"]:

    rule run_bedpost_gpu:
        input:
            diff_dir=rules.copy_inputs_for_bedpost.output.diff_dir,
            dwi=rules.copy_inputs_for_bedpost.output.dwi,
            bval=rules.copy_inputs_for_bedpost.output.bval,
            bvec=rules.copy_inputs_for_bedpost.output.bvec,
            brainmask=rules.copy_inputs_for_bedpost.output.brainmask,
        params:
            container=config["singularity"]["fsl_gpu"],
        output:
            bedpost_dir=directory(
                bids(
                    root=work,
                    desc="eddy",
                    suffix="diffusion.bedpostX",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
        group:
            "subj"
        threads: 32  #needs to be set to avoid multiple gpus from executing
        resources:
            gpus=1,
            mem_mb=16000,
            runtime=360,
        shell:
            #remove the logs to reduce # of files   
            "singularity exec --nv --home $PWD -e {params.container} "
            "bedpostx_gpu {input.diff_dir} && "
            "rm -rf {output.bedpost_dir}/logs "


else:

    rule run_bedpost_cpu:
        input:
            diff_dir=rules.copy_inputs_for_bedpost.output.diff_dir,
            dwi=rules.copy_inputs_for_bedpost.output.dwi,
            bval=rules.copy_inputs_for_bedpost.output.bval,
            bvec=rules.copy_inputs_for_bedpost.output.bvec,
            brainmask=rules.copy_inputs_for_bedpost.output.brainmask,
        params:
            container=config["singularity"]["fsl_cpu"],
            bedpost_script=os.path.join(
                workflow.basedir, "scripts/bedpost/bedpostx-parallel"
            ),
            parallel_script=os.path.join(workflow.basedir, "scripts/parallel"),
            parallel_script_container="/usr/bin/parallel",
        output:
            bedpost_dir=directory(
                bids(
                    root=work,
                    desc="eddy",
                    suffix="diffusion.bedpostX",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
        group:
            "subj"
        threads: 32
        resources:
            mem_mb=16000,
            runtime=360,
        shell:
            "singularity exec --home $PWD "
            "-B {params.parallel_script}:{params.parallel_script_container} "
            "-e {params.container} {params.bedpost_script} "
            "{input.diff_dir} -P {threads} && "
            "rm -rf {output.bedpost_dir}/logs "


rule cp_bedpost_to_results:
    input:
        bedpost_dir=bids(
            root=work,
            desc="eddy",
            suffix="diffusion.bedpostX",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        bedpost_dir=directory(
            bids(
                root=root,
                desc="eddy",
                suffix="diffusion.bedpostX",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            )
        ),
    group:
        "subj"
    shell:
        "cp -Rv {input} {output}"
