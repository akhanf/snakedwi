rule qc:
    input:
        mask_qc=rules.compile_qc_b0_brainmask_manifest.output[0],
        reg_qc=(
            rules.compile_qc_reg_dwi_t1_manifest.output[0]
            if 'T1w' in config.get("output_spaces", ['T1w'])
            else []
        ),
    output:
        os.path.join(qc, "data.json"),
    run:
        qc_data = {
            "mask": json.loads(Path(input["mask_qc"]).read_text()),
        }
        # Only add reg QC if T1w space is selected
        if 'T1w' in config.get("output_spaces", ['T1w']):
            qc_data["reg"] = json.loads(Path(input["reg_qc"]).read_text())

        with open(output[0], "w") as f:
            json.dump(qc_data, f)


_qc_app = os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz")


def _get_tar_contents(tar_path):
    import tarfile

    with tarfile.open(tar_path, "r:gz") as tar:
        return [member.name for member in tar.getmembers() if member.isfile()]


rule unpack_qc_app:
    input:
        os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz"),
    output:
        _get_tar_contents(_qc_app),
    shell:
        "tar -xvzf {input}"
