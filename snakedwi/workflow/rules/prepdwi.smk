
wildcard_constraints:
    shell="[0-9]+",


checkpoint check_subj_dwi_metadata:
    input:
        dwi_jsons=lambda wildcards: expand(
            re.sub(".nii.gz", ".json", input_path["dwi"]),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        workflowopts=directory(
            bids(root=root, datatype="dwi", suffix="workflowopts", **subj_wildcards)
        ),
        metadata=bids(
            root=root, datatype="dwi", suffix="metadata.tsv", **subj_wildcards
        ),
    group:
        "subj"
    script:
        "../scripts/check_subj_dwi_metadata.py"


rule concat_subj_metadata:
    input:
        tsvs=expand(
            bids(root=root, datatype="dwi", suffix="metadata.tsv", **subj_wildcards),
            zip,
            **subj_zip_list
        ),
    output:
        tsv=bids(root=root, suffix="metadata.tsv"),
    container:
        config["singularity"]["python"]
    script:
        "../scripts/concat_tsv.py"


rule create_missing_subj_tsv:
    """creates a tsv file containing subjects that are 
    skipped because they either don't have T1w or don't have dwi data"""
    params:
        missing_subject_zip_list=missing_subj_zip_list,
    output:
        tsv=bids(root=root, suffix="missing.tsv"),
    container:
        config["singularity"]["python"]
    script:
        "../scripts/create_missing_subj_tsv.py"


rule import_dwi:
    input:
        dwi=re.sub(".nii.gz", ".{ext}", input_path["dwi"]),
    output:
        dwi=bids(
            root=work,
            suffix="dwi.{ext,nii.gz|bval|bvec|json}",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    shell:
        "cp {input.dwi} {output.dwi}"


rule dwidenoise:
    input:
        multiext(
            bids(root=work, suffix="dwi", datatype="dwi", **input_wildcards["dwi"]),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    output:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                desc="denoise",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="denoise.log", **input_wildcards["dwi"]),
    group:
        "subj"
    shell:
        "dwidenoise {input[0]} {output[0]} 2> {log} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


def get_degibbs_inputs(wildcards):
    # if input dwi at least 30 dirs, then grab denoised as input
    # else grab without denoising
    import numpy as np

    in_dwi_bval = re.sub(".nii.gz", ".bval", input_path["dwi"].format(**wildcards))
    bvals = np.loadtxt(in_dwi_bval)
    if bvals.size < 30:
        prefix = bids(root=work, suffix="dwi", datatype="dwi", **wildcards)
    else:
        prefix = bids(
            root=work, suffix="dwi", datatype="dwi", desc="denoise", **wildcards
        )
    return multiext(prefix, ".nii.gz", ".bvec", ".bval", ".json")


rule mrdegibbs:
    input:
        get_degibbs_inputs,
    output:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="degibbs.log", **input_wildcards["dwi"]),
    group:
        "subj"
    shell:
        "mrdegibbs {input[0]} {output[0]} 2> {log} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


rule moco_bzeros:
    """ run motion-correction (rigid reg to init volume) on the bzeros """
    input:
        nii_4d=bids(
            root=work,
            suffix="b0s.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    output:
        affine_dir=directory(
            bids(
                root=work,
                suffix="transforms",
                desc="moco",
                datatype="dwi",
                **subj_wildcards
            )
        ),
        nii_4d=bids(
            root=work,
            suffix="b0s.nii.gz",
            desc="moco",
            datatype="dwi",
            **subj_wildcards
        ),
        nii_avg3d=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="moco",
            datatype="dwi",
            **subj_wildcards
        ),
    threads: 8  #doesn't need to be more than the number of bzeros 
    resources:
        mem_mb=32000,
    shadow:
        "minimal"
    container:
        config["singularity"]["prepdwi"]  #-- this rule needs niftyreg, c3d and mrtrix
    group:
        "subj"
    shell:
        "c4d {input.nii_4d} -slice w 0:-1 -oo dwi_%03d.nii && "
        "parallel --eta --jobs {threads} "
        "reg_aladin -flo dwi_{{1}}.nii  -ref dwi_000.nii -res warped_{{1}}.nii -aff affine_xfm_ras_{{1}}.txt "
        " ::: `ls dwi_???.nii | tail -n +2 | grep -Po '(?<=dwi_)[0-9]+'` && "
        " mkdir -p {output.affine_dir} && cp affine_xfm_ras_*.txt {output.affine_dir} && "
        " echo -e '1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1' > {output.affine_dir}/affine_xfm_ras_000.txt && "
        " mrcat dwi_000.nii warped_*.nii {output.nii_4d} && "
        " mrmath {output.nii_4d} mean {output.nii_avg3d} -axis 3"


# now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
        json=bids(
            root=work,
            suffix="dwi.json",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
    output:
        phenc_txt=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/get_phase_encode_txt.py"


rule concat_phase_encode_txt:
    input:
        phenc_txts=lambda wildcards: expand(
            bids(
                root=work,
                suffix="phenc.txt",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        phenc_concat=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cat {input} > {output}"


# function for either concatenating (if multiple inputs) or copying
def get_concat_or_cp_cmd(wildcards, input, output):
    if len(input) > 1:
        cmd = f"mrcat {input} {output}"
    elif len(input) == 1:
        cmd = f"cp {input} {output}"
    else:
        # no inputs
        cmd = None
    return cmd


rule concat_bzeros:
    input:
        bzero_niis=lambda wildcards: expand(
            bids(
                root=work,
                suffix="b0.nii.gz",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        bzero_concat=bids(
            root=work,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="concat_bzeros.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


rule concat_degibbs_dwi:
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
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        dwi_concat=bids(
            root=work,
            suffix="dwi.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="concat_degibbs_dwi.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


rule concat_runs_bvec:
    input:
        lambda wildcards: expand(
            bids(
                root=work,
                suffix="dwi.bvec",
                desc="{{desc}}",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        bids(
            root=work,
            suffix="dwi.bvec",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/concat_bv.py"


rule concat_runs_bval:
    input:
        lambda wildcards: expand(
            bids(
                root=work,
                suffix="dwi.bval",
                desc="{{desc}}",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        bids(
            root=work,
            suffix="dwi.bval",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/concat_bv.py"


# combines json files from multiple scans -- for now as a hack just copying first json over..
rule concat_runs_json:
    input:
        lambda wildcards: expand(
            bids(
                root=work,
                suffix="dwi.json",
                desc="{{desc}}",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    output:
        bids(
            root=work,
            suffix="dwi.json",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input[0]} {output}"


#    script: '../scripts/concat_json.py'


rule get_shells_from_bvals:
    input:
        "{dwi_prefix}.bval",
    output:
        "{dwi_prefix}.shells.json",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/get_shells_from_bvals.py"


# writes 4d file
rule get_shell_avgs:
    input:
        dwi="{dwi_prefix}.nii.gz",
        shells="{dwi_prefix}.shells.json",
    output:
        avgshells="{dwi_prefix}.avgshells.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/get_shell_avgs.py"


# this gets a particular shell (can use to get b0)
rule get_shell_avg:
    input:
        dwi="{dwi_prefix}_dwi.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        avgshell="{dwi_prefix}_b{shell}.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/get_shell_avg.py"


# this gets vols from a a particular shell (can use to get b0s)
rule get_shell_vols:
    input:
        dwi="{dwi_prefix}_dwi.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        avgshell="{dwi_prefix}_b{shell}s.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/get_shell_vols.py"


def get_b0_init_mask():
    return bids(
        root=work,
        suffix="mask.nii.gz",
        desc="brain",
        method="synthstrip",
        from_="mocob0",
        datatype="dwi",
        **subj_wildcards,
    )


def get_b0_mask():
    if config["masking_method"] == "b0_BET":
        method = "bet_from-b0"
    elif config["masking_method"] == "b0_SyN":
        method = f"b0SyN_from-{config['template']}"
    elif config["masking_method"] == "b0_synthstrip":
        method = "synthstrip_from-dwirefb0"

    # then get bids name of file
    return bids(
        root=work,
        suffix="mask.nii.gz",
        desc="brain",
        method=method,
        datatype="dwi",
        **subj_wildcards,
    )


# generate qc snapshot for brain  mask
rule qc_b0_brainmask:
    input:
        img=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="dwiref",
            datatype="dwi",
            **subj_wildcards
        ),
        seg=get_b0_mask(),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                suffix="mask.png",
                desc="brain",
                **subj_wildcards
            ),
            caption="../report/brainmask_dwi.rst",
            category="Brainmask",
        ),
        html=bids(
            root=root,
            datatype="qc",
            suffix="mask.html",
            desc="brain",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/vis_qc_dseg.py"
