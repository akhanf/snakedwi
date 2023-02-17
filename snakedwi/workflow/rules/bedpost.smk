
rule copy_inputs_for_bedpost:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="eddy",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bval=bids(
            root=root,
            suffix="dwi.bval",
            desc="eddy",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bvec=bids(
            root=root,
            suffix="dwi.bvec",
            desc="eddy",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        brainmask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
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
            diff_dir=bids(
                root=work,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
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
        params:
            container=config["singularity"]["bedpost_gpu"],
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
        threads: 32  #this needs to be set in order to avoid multiple gpus from executing
        resources:
            gpus=1,
            mem_mb=16000,
            time=360,
        shell:
            #remove the logs to reduce # of files  
            # remove the input dir (copy of files) 
            "singularity exec --nv --home $PWD -e {params.container} bedpostx_gpu {input.diff_dir} && "
            "rm -rf {output.bedpost_dir}/logs "


else:

    rule run_bedpost_cpu:
        input:
            diff_dir=bids(
                root=work,
                desc="eddy",
                suffix="diffusion",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
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
        params:
            container=config["singularity"]["fsl_cpu"],
            bedpost_script=os.path.join(workflow.basedir, "scripts/bedpostx-parallel"),
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
            time=360,
        shell:
            "singularity exec --home $PWD -B {params.parallel_script}:{params.parallel_script_container} -e "
            " {params.container} {params.bedpost_script} {input.diff_dir} -P {threads} && "
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
