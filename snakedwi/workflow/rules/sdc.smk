rule run_topup:
    input:
        bzero_concat=bids(
            root=work,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
        phenc_concat=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    params:
        out_prefix=bids(root=work, suffix="topup", datatype="dwi", **subj_wildcards),
        config="b02b0.cnf",  #this config sets the multi-res schedule and other params..
    output:
        bzero_corrected=bids(
            root=work,
            suffix="concatb0.nii.gz",
            desc="topup",
            datatype="dwi",
            **subj_wildcards
        ),
        fieldmap=bids(
            root=work,
            suffix="fmap.nii.gz",
            desc="topup",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_fieldcoef=bids(
            root=work,
            suffix="topup_fieldcoef.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_movpar=bids(
            root=work, suffix="topup_movpar.txt", datatype="dwi", **subj_wildcards
        ),
    container:
        config["singularity"]["fsl"]
    log:
        bids(root="logs", suffix="topup.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "topup --imain={input.bzero_concat} --datain={input.phenc_concat} --config={params.config}"
        " --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v 2> {log}"


# this is for equal positive and negative blipped data - method=lsr --unused currently (jac method can be applied more generally)
rule apply_topup_lsr:
    input:
        dwi_niis=lambda wildcards: expand(
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
        phenc_concat=bids(
            root=work,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_fieldcoef=bids(
            root=work,
            suffix="topup_fieldcoef.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_movpar=bids(
            root=work, suffix="topup_movpar.txt", datatype="dwi", **subj_wildcards
        ),
    params:
        #create comma-seperated list of dwi nii
        imain=lambda wildcards, input: ",".join(input.dwi_niis),
        # create comma-sep list of indices 1-N
        inindex=lambda wildcards, input: ",".join(
            [str(i) for i in range(1, len(input.dwi_niis) + 1)]
        ),
        topup_prefix=bids(root=work, suffix="topup", datatype="dwi", **subj_wildcards),
        out_prefix="dwi_topup",
    output:
        dwi_topup=bids(
            root=work,
            suffix="dwi.nii.gz",
            desc="topup",
            method="lsr",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["fsl"]
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        "applytopup --verbose --datain={input.phenc_concat} --imain={params.imain} --inindex={params.inindex} "
        " -t {params.topup_prefix} -o {params.out_prefix} && "
        " fslmaths {params.out_prefix}.nii.gz {output.dwi_topup}"


def get_applytopup_inindex(wildcards):

    # need to get index of the scan from the subject scans
    # so first filter the ziplist to get the subject(+session)
    subj_filter = {"subject": wildcards.subject}
    if "session" in subj_wildcards.keys():
        subj_filter["session"] = wildcards.session

    zip_list_subj = filter_list(input_zip_lists["dwi"], subj_filter)

    # now filter the subj ziplist using all wildcards to get the index of the scan
    indices = filter_list(zip_list_subj, wildcards, return_indices_only=True)
    return indices[0] + 1  # get first index, but adjust to start at 1 instead of 0


rule apply_topup_jac:
    input:
        nii=bids(
            root=work,
            suffix="dwi.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
        phenc_scan=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
        phenc_concat=bids(
            root=work,
            suffix="phenc.txt",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_fieldcoef=bids(
            root=work,
            suffix="topup_fieldcoef.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
        topup_movpar=bids(
            root=work, suffix="topup_movpar.txt", datatype="dwi", **subj_wildcards
        ),
    params:
        inindex=get_applytopup_inindex,
        topup_prefix=bids(root=work, suffix="topup", datatype="dwi", **subj_wildcards),
    output:
        nii=bids(
            root=work,
            suffix="dwi.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    container:
        config["singularity"]["fsl"]
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        " applytopup --verbose --datain={input.phenc_concat} --imain={input.nii} --inindex={params.inindex} "
        " -t {params.topup_prefix} -o dwi_topup --method=jac && mv dwi_topup.nii.gz {output.nii}"


# topup-corrected data is only used for brainmasking..
# here, use the jac method by default (later can decide if lsr approach can be used based on headers)
# with jac approach, the jac images need to be concatenated, then avgshell extracted

"""
rule cp_sidecars_topup_lsr:
    #TODO: BEST WAY TO TO EXEMPLAR DWI? 
    input: multiext(bids(root='work',suffix='dwi',desc='degibbs',datatype='dwi',**config['subj_wildcards'],**dwi_exemplar_dict),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='work',suffix='dwi',desc='topup',method='lsr',datatype='dwi',**config['subj_wildcards']),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')
"""


rule cp_sidecars_topup_jac:
    input:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                desc="degibbs",
                datatype="dwi",
                **subj_wildcards
            ),
            ".bvec",
            ".bval",
            ".json",
        ),
    output:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                desc="topup",
                method="jac",
                datatype="dwi",
                **subj_wildcards
            ),
            ".bvec",
            ".bval",
            ".json",
        ),
    group:
        "subj"
    run:
        for in_file, out_file in zip(input, output):
            shell("cp -v {in_file} {out_file}")


rule concat_dwi_topup_jac:
    input:
        dwi_niis=lambda wildcards: expand(
            bids(
                root=work,
                suffix="dwi.nii.gz",
                desc="topup",
                method="jac",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        dwi_concat=bids(
            root=work,
            suffix="dwi.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    group:
        "subj"
    shell:
        "mrcat {input} {output}"


rule syn_sdc:
    input:
        in_epis=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
        in_anat=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
        json_file=lambda wildcards: expand(
            bids(
                root=work, suffix="dwi.json", datatype="dwi", **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        )[0],
        mask_anat=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="brain",
            suffix="mask.nii.gz"
        ),
        epi_mask=get_b0_init_mask(),
        std2anat_xfm=bids(
            root=work,
            datatype="transforms",
            from_=config["template"],
            to="subject",
            desc="affine",
            type_="itk",
            **subj_wildcards,
            suffix="xfm.txt"
        ),
    output:
        unwarped=bids(
            root=work,
            datatype="dwi",
            suffix="b0.nii.gz",
            desc="unwarped",
            method="synsdc",
            **subj_wildcards
        ),
        unwarped_mask=bids(
            root=work,
            datatype="dwi",
            suffix="mask.nii.gz",
            desc="brain",
            method="synsdc",
            **subj_wildcards
        ),
        xfm=bids(
            root=work,
            datatype="dwi",
            suffix="xfm.nii.gz",
            desc="itk",
            method="synsdc",
            **subj_wildcards
        ),
        fmap=bids(
            root=work,
            datatype="dwi",
            suffix="fmap.nii.gz",
            desc="b0",
            method="synsdc",
            **subj_wildcards
        ),
    threads: 8
    group:
        "subj"
    shadow:
        "minimal"
    container:
        config["singularity"]["sdcflows"]
    script:
        "../scripts/sdcflows_syn.py"


def get_dwi_ref(wildcards):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([method],) = glob_wildcards(os.path.join(checkpoint_output, "sdc-{method}"))

    if method == "topup":
        return bids(
            root=work,
            suffix="b0.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **subj_wildcards
        )
    else:
        return bids(
            root=work,
            datatype="dwi",
            suffix="b0.nii.gz",
            desc="unwarped",
            method="synsdc",
            **subj_wildcards
        )


rule cp_dwi_ref:
    input:
        get_dwi_ref,
    output:
        bids(
            root=work,
            suffix="b0.nii.gz",
            desc="dwiref",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"
