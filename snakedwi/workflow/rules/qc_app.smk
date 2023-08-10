rule qc:
    input:
        mask_qc=rules.compile_qc_b0_brainmask_manifest.output[0],
        reg_qc=rules.compile_qc_reg_dwi_t1_manifest.output[0],
    output:
        os.path.join(qc, "data.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "mask": json.loads(Path(input["mask_qc"]).read_text()),
                    "reg": json.loads(Path(input["reg_qc"]).read_text()),
                },
                f,
            )


_qc_app = os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz")


def _get_tar_contents(file):
    try:
        return [
            p
            for p in sp.check_output(["tar", "-tf", _qc_app]).decode().splitlines()
            if p[-1] != "/"
        ]
    except sp.CalledProcessError as err:
        raise Exception("Unable to find qc-app.tar.gz...") from err


rule unpack_qc_app:
    input:
        os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz"),
    output:
        _get_tar_contents(_qc_app),
    shell:
        "tar -xvzf {input}"
