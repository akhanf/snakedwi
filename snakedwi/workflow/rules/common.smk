def get_eddy_quad_all():
    if config["eddy_no_quad"]:
        return {}
    else:
        return {
            "eddy_qc": expand(
                bids(root=root, datatype="qc", suffix="eddyqc", **subj_wildcards),
                zip,
                **subj_zip_list,
            )
        }


def get_bedpost_all():
    # Only run bedpost if T1w space is in output_spaces
    if config["no_bedpost"] or 'T1w' not in config.get("output_spaces", ['T1w']):
        return {}
    else:
        return {
            "bedpost": expand(
                bids(
                    root=root,
                    datatype="dwi",
                    suffix="diffusion.bedpostX",
                    desc="eddy",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    **subj_wildcards,
                ),
                zip,
                **subj_zip_list,
            )
        }


def get_dtifit_t1w_all():
    """Get dtifit outputs in T1w space"""
    if 'T1w' not in config.get("output_spaces", ['T1w']):
        return {}
    else:
        return {
            "dtifit_t1w": expand(
                bids(
                    root=root,
                    datatype="dwi",
                    suffix="dtifit",
                    desc="eddy",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    **subj_wildcards,
                ),
                zip,
                **subj_zip_list,
            )
        }


def get_dtifit_dwi_all():
    """Get dtifit outputs in dwi (native) space"""
    if 'dwi' not in config.get("output_spaces", ['T1w']):
        return {}
    else:
        return {
            "dtifit_dwi": expand(
                bids(
                    root=root,
                    datatype="dwi",
                    suffix="dtifit",
                    desc="eddy",
                    **subj_wildcards,
                ),
                zip,
                **subj_zip_list,
            )
        }
