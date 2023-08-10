def _get_mocorrected_b0s(wcards):
    b0s = get_dwi_indices(
        expand(
            bids(
                root=work,
                suffix="b0.nii.gz",
                datatype="dwi",
                desc="moco",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wcards)
        ),
        wcards,
    )
    if len(b0s) == 1:
        return b0s
    return rules.moco_bzeros_3d.output["nii_4d"]


rule run_topup:
    input:
        bzero_concat=rules.concat_bzeros.output.bzero_concat,
        phenc_concat=rules.concat_phase_encode_txt.output.phenc_concat,
    params:
        out_prefix=bids(
            root=work,
            suffix="topup",
            datatype="dwi",
            **subj_wildcards,
        ),
        config="b02b0.cnf",  #sets the multi-res schedule and other params..
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
            root=work,
            suffix="topup_movpar.txt",
            datatype="dwi",
            **subj_wildcards,
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


def get_applytopup_inindex(wildcards):
    index = get_index_of_dwi_scan(wildcards)
    return index + 1  # adjust to start at 1 instead of 0


rule apply_topup_jac:
    input:
        nii=rules.moco_scan_bzeros_4d.output.nii_avg3d,
        phenc_scan=rules.get_phase_encode_txt.output.phenc_txt,
        phenc_concat=rules.concat_phase_encode_txt.output.phenc_concat,
        topup_fieldcoef=rules.run_topup.output.topup_fieldcoef,
        topup_movpar=rules.run_topup.output.topup_movpar,
    params:
        inindex=get_applytopup_inindex,
        topup_prefix=bids(
            root=work,
            suffix="topup",
            datatype="dwi",
            **subj_wildcards,
        ),
    output:
        nii=bids(
            root=work,
            suffix="b0.nii.gz",
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


ruleorder: avg_b0s_topup_jac > get_shell_avg


rule avg_b0s_topup_jac:
    input:
        b0_niis=lambda wildcards: get_dwi_indices(
            expand(
                rules.apply_topup_jac.output.nii,
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        b0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="topup",
            method="jac",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        "mrcat {input} b0_concat.nii && "
        "mrmath b0_concat.nii mean {output} -axis 3"
