checkpoint check_subj_dwi_metadata:
    input:
        dwi_jsons=lambda wildcards: expand(
            re.sub(".nii.gz", ".json", input_path["dwi"]),
            zip,
            **filter_list(input_zip_lists["dwi"], wildcards)
        ),
    params:
        index_col_value=bids(
            include_subject_dir=False,
            include_session_dir=False,
            **subj_wildcards,
        ),
        index_col_name="subj",
    output:
        workflowopts=directory(
            bids(
                root=root,
                datatype="dwi",
                suffix="workflowopts",
                **subj_wildcards,
            )
        ),
        metadata=bids(
            root=root,
            datatype="dwi",
            suffix="metadata.tsv",
            **subj_wildcards,
        ),
    script:
        "../scripts/metadata/check_subj_dwi_metadata.py"


rule concat_subj_metadata:
    input:
        tsvs=expand(
            rules.check_subj_dwi_metadata.output.metadata,
            zip,
            **subj_zip_list,
        ),
    output:
        tsv=bids(root=root, suffix="metadata.tsv"),
    container:
        config["singularity"]["python"]
    shell:
        "../scripts/metadata/concat_subj_metadata.py"


rule create_missing_subj_tsv:
    """Create tsv file of skipped subjects due to missing T1w / dwi data"""
    params:
        missing_subject_zip_list=missing_subj_zip_list,
    output:
        tsv=bids(root=root, suffix="missing.tsv"),
    container:
        config["singularity"]["python"]
    script:
        "../scripts/metadata/create_missing_subj_tsv.py"


def get_dwi_indices(all_dwi, wildcards):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([indices_str],) = glob_wildcards(
        os.path.join(checkpoint_output, "indices-{indices}")
    )
    indices = indices_str.split(",")
    return [all_dwi[int(i)] for i in indices]


def get_dwi_indices_only(wildcards):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([indices_str],) = glob_wildcards(
        os.path.join(checkpoint_output, "indices-{indices}")
    )
    indices = indices_str.split(",")
    return [int(i) for i in indices]


def get_sdc_method(wildcards):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    # this gets the sdc method for this subject(session)
    ([method],) = glob_wildcards(os.path.join(checkpoint_output, "sdc-{method}"))
    return method


def get_dwi_num_scans(wildcards):
    # this gets the number of DWI scans for this subject(session)
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([indices_str],) = glob_wildcards(
        os.path.join(checkpoint_output, "indices-{indices}")
    )
    indices = indices_str.split(",")
    return len(indices)


def get_pe_axis(wildcards):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([peaxis],) = glob_wildcards(os.path.join(checkpoint_output, "PEaxis-{axis}"))
    return peaxis


def get_enable_s2v(wildcards):
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([s2v_is_enabled],) = glob_wildcards(
        os.path.join(checkpoint_output, "eddys2v-{isenabled}")
    )
    return s2v_is_enabled


def get_index_of_dwi_scan(wildcards):
    """given wildcards into a specific dwi acquisition, this returns the 0-based
    index of that scan in the list of dwi scans in use for this acquisition"""
    checkpoint_output = checkpoints.check_subj_dwi_metadata.get(**wildcards).output[0]
    ([indices_str],) = glob_wildcards(
        os.path.join(checkpoint_output, "indices-{indices}")
    )
    dwi_indices = [int(i) for i in indices_str.split(",")]

    # first filter the ziplist to get the subject(+session)
    subj_filter = {"subject": wildcards.subject}
    if "session" in subj_wildcards.keys():
        subj_filter["session"] = wildcards.session

    zip_list_subj = filter_list(input_zip_lists["dwi"], subj_filter)

    # filter the zip list using the indices
    for zkey, zlist in zip_list_subj.items():
        zip_list_subj[zkey] = [zlist[i] for i in dwi_indices]

    # now filter the subj ziplist using all wildcards to get the index of the scan
    indices = filter_list(zip_list_subj, wildcards, return_indices_only=True)
    return indices[0]
