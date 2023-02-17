def get_eddy_quad_all():
    if config["eddy_no_quad"]:
        return {}
    else:
        return {
            "eddy_qc": expand(
                bids(root=root, datatype="qc", suffix="eddyqc", **subj_wildcards),
                zip,
                **input_zip_lists["dwi"]
            )
        }


def get_bedpost_all():
    if config["no_bedpost"]:
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
                    **subj_wildcards
                ),
                zip,
                **input_zip_lists["dwi"]
            )
        }
