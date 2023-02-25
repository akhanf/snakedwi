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
        dwi_niis=lambda wildcards: get_dwi_indices(
            expand(
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
            wildcards,
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
        json_file=lambda wildcards: get_dwi_indices(
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
        "python /SynthSR/scripts/predict_command_line.py  --cpu --threads {threads} {input} {output}"


rule reslice_synthSR_b0:
    input:
        ref=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
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
            **subj_wildcards,
            desc="preproc",
            suffix="T1wSynthSR.nii.gz"
        ),
        ref=bids(
            root=work,
            suffix="b0SynthSRresliced.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
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
        bids(root="logs", suffix="rigid_reg_t1_to_b0_synthsr.txt", **subj_wildcards),
    group:
        "subj"
    shell:
        "reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.xfm_ras} > {log}"


def get_restrict_deformation(wildcards, input, output):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([peaxis],) = glob_wildcards(os.path.join(checkpoint_output, "PEaxis-{axis}"))
    if peaxis == "i":
        return "1x0x0"
    elif peaxis == "j":
        return "0x1x0"
    elif peaxis == "k":
        return "0x0x1"  # unlikely..


rule reg_b0_to_t1_synthsr:
    input:
        rules.check_subj_dwi_metadata.output,
        t1synth=bids(
            root=work,
            suffix="T1wSynthSRreg.nii.gz",
            space="rigidb0",
            datatype="dwi",
            **subj_wildcards
        ),
        b0synth=bids(
            root=work,
            suffix="b0SynthSRresliced.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
    params:
        general_opts="-d 3 -v",
        metric=lambda wildcards, input: f"MeanSquares[{input.t1synth},{input.b0synth}]",
        transform="SyN[0.1]",
        convergence="100x50x20",
        smoothing_sigmas="2x1x0vox",
        shrink_factors="4x2x1",
        restrict_deformation=get_restrict_deformation,
    # "0x1x0",  #should be set to the phase encode dir - can read json for this..
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
        "antsRegistration {params.general_opts} --metric {params.metric} --convergence {params.convergence} "
        "  --smoothing-sigmas {params.smoothing_sigmas} --shrink-factors {params.shrink_factors}  "
        "  --transform {params.transform} "
        " --output [ants_,{output.unwarped}]  --restrict-deformation {params.restrict_deformation} > {log} && "
        " mv ants_0Warp.nii.gz {output.fwd_xfm} "


rule displacement_field_to_fmap:
    input:
        disp_field=bids(
            root=root,
            suffix="xfm.nii.gz",
            from_="b0SynthSR",
            to="T1wSynthSR",
            type_="itk",
            desc="SyN",
            datatype="transforms",
            **subj_wildcards
        ),
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
        "../scripts/displacement_field_to_fmap.py"


rule apply_unwarp_synthsr:
    input:
        b0=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="moco",
            **subj_wildcards
        ),
        xfm=bids(
            root=root,
            suffix="xfm.nii.gz",
            from_="b0SynthSR",
            to="T1wSynthSR",
            type_="itk",
            desc="SyN",
            datatype="transforms",
            **subj_wildcards
        ),
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
        " -r {input.b0} -t {input.xfm} -e 0 -n NearestNeighbor"


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
    elif method == "synthsr":
        return bids(
            root=work,
            datatype="dwi",
            suffix="b0.nii.gz",
            desc="unwarped",
            method="synthSRsdc",
            **subj_wildcards
        )

    elif method == "syn":
        return bids(
            root=work,
            datatype="dwi",
            suffix="b0.nii.gz",
            desc="unwarped",
            method="synsdc",
            **subj_wildcards
        )
    else:
        return (
            bids(
                root=work,
                suffix="b0.nii.gz",
                datatype="dwi",
                desc="moco",
                **subj_wildcards
            ),
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
