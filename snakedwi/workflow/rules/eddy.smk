rule get_eddy_index_txt:
    input:
        dwi_niis=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.nii.gz",
                    desc="degibbs",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        eddy_index_txt=bids(
            root=work,
            suffix="dwi.eddy_index.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/get_eddy_index_txt.py"


# if the --slspec_txt option is not used, use the SliceTiming json, otherwise just get from the file supplied at command-line
if not config["slspec_txt"]:

    rule get_slspec_txt:
        input:
            dwi_jsons=lambda wildcards: get_dwi_indices(
                expand(
                    bids(
                        root=work,
                        suffix="dwi.json",
                        desc="degibbs",
                        datatype="dwi",
                        **input_wildcards["dwi"]
                    ),
                    zip,
                    **filter_list(input_zip_lists["dwi"], wildcards)
                ),
                wildcards,
            ),
        output:
            eddy_slspec_txt=bids(
                root=work,
                suffix="dwi.eddy_slspec.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
        group:
            "subj"
        container:
            config["singularity"]["python"]
        script:
            "../scripts/get_slspec_txt.py"


else:

    rule get_slspec_txt:
        input:
            config["slspec_txt"],
        output:
            eddy_slspec_txt=bids(
                root=work,
                suffix="dwi.eddy_slspec.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
        group:
            "subj"
        shell:
            "cp {input} {output}"


def get_eddy_topup_fmap_input(wildcards):
    # this gets the number of DWI scans for this subject(session)
    filtered = filter_list(input_zip_lists["dwi"], wildcards)
    num_scans = len(filtered["subject"])

    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([method],) = glob_wildcards(os.path.join(checkpoint_output, "sdc-{method}"))

    if num_scans > 1 and method == "topup":
        return {
            "topup_fieldcoef": bids(
                root=work,
                suffix=f"topup_fieldcoef.nii.gz",
                datatype="dwi",
                **subj_wildcards,
            ).format(**wildcards),
            "topup_movpar": bids(
                root=work,
                suffix=f"topup_movpar.txt",
                datatype="dwi",
                **subj_wildcards,
            ).format(**wildcards),
        }
    elif method == "synthsr":
        return {
            "fmap": bids(
                root=work,
                datatype="dwi",
                suffix="fmap.nii.gz",
                desc="b0",
                method="synthSRsdc",
                **subj_wildcards,
            ).format(**wildcards)
        }
    elif method == "syn":
        return {
            "fmap": bids(
                root=work,
                datatype="dwi",
                suffix="fmap.nii.gz",
                desc="b0",
                method="synsdc",
                **subj_wildcards,
            ).format(**wildcards)
        }
    else:
        return {}


def get_eddy_topup_fmap_opt(wildcards, input):

    # this gets the number of DWI scans for this subject(session)
    filtered = filter_list(input_zip_lists["dwi"], wildcards)
    num_scans = len(filtered["subject"])

    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([method],) = glob_wildcards(os.path.join(checkpoint_output, "sdc-{method}"))

    if num_scans > 1 and method == "topup":
        topup_prefix = bids(
            root=work, suffix="topup", datatype="dwi", **subj_wildcards
        ).format(**wildcards)
        return f"--topup={topup_prefix}"
    elif method == "synthsr":
        fmap_prefix = bids(
            root=work,
            datatype="dwi",
            suffix="fmap",
            desc="b0",
            method="synthSRsdc",
            **subj_wildcards,
        ).format(**wildcards)
        return f"--field={fmap_prefix}"

    elif method == "syn":
        fmap_prefix = bids(
            root=work,
            datatype="dwi",
            suffix="fmap",
            desc="b0",
            method="synsdc",
            **subj_wildcards,
        ).format(**wildcards)
        return f"--field={fmap_prefix}"
    else:
        return ""


def get_eddy_s2v_opts(wildcards, input):

    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([s2v_is_enabled],) = glob_wildcards(
        os.path.join(checkpoint_output, "eddys2v-{isenabled}")
    )

    options = []
    if s2v_is_enabled == "yes":
        options += [
            f"--{key}={value}"
            for (key, value) in config["eddy"]["with_s2v"].items()
            if value is not None
        ]
    else:
        options += [
            f"--{key}={value}"
            for (key, value) in config["eddy"]["without_s2v"].items()
            if value is not None
        ]

    return " ".join(options)


def get_eddy_slspec_input(wildcards):

    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([s2v_is_enabled],) = glob_wildcards(
        os.path.join(checkpoint_output, "eddys2v-{isenabled}")
    )

    if s2v_is_enabled == "yes":
        return {
            "eddy_slspec_txt": bids(
                root=work,
                suffix="dwi.eddy_slspec.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards,
            )
        }

    else:
        return {}


def get_eddy_slspec_opt(wildcards, input):

    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([s2v_is_enabled],) = glob_wildcards(
        os.path.join(checkpoint_output, "eddys2v-{isenabled}")
    )

    if s2v_is_enabled == "yes":
        return f"--slspec={input.eddy_slspec_txt}"
    else:
        return ""


if config["use_eddy_gpu"]:

    rule run_eddy_gpu:
        input:
            unpack(get_eddy_slspec_input),
            unpack(get_eddy_topup_fmap_input),
            dwi_concat=bids(
                root=work,
                suffix="dwi.nii.gz",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            phenc_concat=bids(
                root=work,
                suffix="phenc.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            eddy_index_txt=bids(
                root=work,
                suffix="dwi.eddy_index.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            brainmask=get_b0_mask(),
            bvals=bids(
                root=work,
                suffix="dwi.bval",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            bvecs=bids(
                root=work,
                suffix="dwi.bvec",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
        params:
            #set eddy output prefix to 'dwi' inside the output folder
            out_prefix=lambda wildcards, output: os.path.join(output.out_folder, "dwi"),
            flags=" ".join(
                [
                    f"--{key}"
                    for (key, value) in config["eddy"]["flags"].items()
                    if value == True
                ]
            ),
            container=config["singularity"]["eddy_gpu"],
            topup_opt=get_eddy_topup_fmap_opt,
            s2v_opts=get_eddy_s2v_opts,
            slspec_opt=get_eddy_slspec_opt,
        output:
            #eddy creates many files, so write them to a eddy subfolder instead
            out_folder=directory(
                bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards)
            ),
            dwi=os.path.join(
                bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
                "dwi.nii.gz",
            ),
            bvec=os.path.join(
                bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
                "dwi.eddy_rotated_bvecs",
            ),
        threads: 16  #this needs to be set in order to avoid multiple gpus from executing
        resources:
            gpus=1,
            time=360,  #6 hours (this is a conservative estimate, may be shorter)
            mem_mb=32000,
        group:
            "subj"
        shell:
            "singularity exec --nv --home $PWD -e {params.container} eddy_cuda9.1 "
            " --imain={input.dwi_concat} --mask={input.brainmask} "
            " --acqp={input.phenc_concat} --index={input.eddy_index_txt} "
            " --bvecs={input.bvecs} --bvals={input.bvals} "
            " --out={params.out_prefix} "
            " {params.s2v_opts} "
            " {params.slspec_opt} "
            " {params.topup_opt} "
            " {params.flags}"


else:

    rule run_eddy_cpu:
        input:
            unpack(get_eddy_slspec_input),
            unpack(get_eddy_topup_fmap_input),
            dwi_concat=bids(
                root=work,
                suffix="dwi.nii.gz",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            phenc_concat=bids(
                root=work,
                suffix="phenc.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            eddy_index_txt=bids(
                root=work,
                suffix="dwi.eddy_index.txt",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            brainmask=get_b0_mask(),
            bvals=bids(
                root=work,
                suffix="dwi.bval",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            bvecs=bids(
                root=work,
                suffix="dwi.bvec",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
        params:
            #set eddy output prefix to 'dwi' inside the output folder
            out_prefix=lambda wildcards, output: os.path.join(output.out_folder, "dwi"),
            flags=" ".join(
                [
                    f"--{key}"
                    for (key, value) in config["eddy"]["flags"].items()
                    if value == True
                ]
            ),
            topup_opt=get_eddy_topup_fmap_opt,
            s2v_opts=get_eddy_s2v_opts,
            slspec_opt=get_eddy_slspec_opt,
        output:
            #eddy creates many files, so write them to a eddy subfolder instead
            out_folder=directory(
                bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards)
            ),
            dwi=os.path.join(
                bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
                "dwi.nii.gz",
            ),
            bvec=os.path.join(
                bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
                "dwi.eddy_rotated_bvecs",
            ),
        threads: 16  #this needs to be set in order to avoid multiple gpus from executing
        resources:
            time=360,  #6 hours (this is a conservative estimate, may be shorter)
            mem_mb=32000,
        log:
            bids(root="logs", suffix="run_eddy.log", **subj_wildcards),
        container:
            config["singularity"]["fsl"]
        group:
            "subj"
        shell:
            "eddy_openmp "
            " --imain={input.dwi_concat} --mask={input.brainmask} "
            " --acqp={input.phenc_concat} --index={input.eddy_index_txt} "
            " --bvecs={input.bvecs} --bvals={input.bvals} "
            " --out={params.out_prefix} "
            " {params.s2v_opts} "
            " {params.slspec_opt} "
            " {params.topup_opt} "
            " {params.flags}  &> {log}"


rule cp_eddy_outputs:
    input:
        #get nii.gz, bvec, and bval from eddy output
        dwi=os.path.join(
            bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
            "dwi.nii.gz",
        ),
        bvec=os.path.join(
            bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
            "dwi.eddy_rotated_bvecs",
        ),
        bval=bids(
            root=work,
            suffix="dwi.bval",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        mask=get_b0_mask(),
    output:
        multiext(
            bids(
                root=root, suffix="dwi", desc="eddy", datatype="dwi", **subj_wildcards
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
        ),
        bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    run:
        for in_file, out_file in zip(input, output):
            shell("cp -v {in_file} {out_file}")




rule eddy_quad:
    input:
        unpack(get_eddy_slspec_input),
        phenc_concat=bids(
            root=work,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        eddy_index_txt=bids(
            root=work,
            suffix="dwi.eddy_index.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        brainmask=get_b0_mask(),
        bvals=bids(
            root=work,
            suffix="dwi.bval",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        bvecs=bids(
            root=work,
            suffix="dwi.bvec",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        eddy_dir=bids(root=work, suffix="eddy", datatype="dwi", **subj_wildcards),
    params:
        eddy_prefix=lambda wildcards, input: os.path.join(input.eddy_dir, "dwi"),
        slspec_opt=get_eddy_slspec_opt,
    output:
        out_dir=directory(
            bids(root=root, suffix="eddyqc", datatype="qc", **subj_wildcards)
        ),
        eddy_qc_pdf=bids(
            root=root, suffix="eddyqc/qc.pdf", datatype="qc", **subj_wildcards
        ),
    container:
        config["singularity"]["fsl"]
    group:
        "subj"
    shell:
        "rmdir {output.out_dir} && "
        "eddy_quad {params.eddy_prefix} -idx {input.eddy_index_txt} -par {input.phenc_concat} "
        " -m {input.brainmask} -b {input.bvals} -g {input.bvecs} -o {output.out_dir} "
        " {params.slspec_opt} -v"
