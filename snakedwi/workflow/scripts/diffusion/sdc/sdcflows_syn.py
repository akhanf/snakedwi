import json
import tempfile
from glob import glob
from pathlib import Path
from shutil import copyfile

import nibabel as nib
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
from nipype import Workflow
from sdcflows.workflows.fit.syn import (
    init_syn_preprocessing_wf,
    init_syn_sdc_wf,
)


def sdcflows_syn(
    in_epis: str,
    in_anat: str,
    mask_anat: str,
    std2anat_xfm: str,
    epi_mask: str,
    json_file: str,
    smk_out: dict,
    threads: int,
) -> None:
    with tempfile.TemporaryDirectory() as tmpdirname:

        workflow = Workflow(name="nipype_wf")

        syn_preprocessing_wf = init_syn_preprocessing_wf(
            omp_nthreads=threads,
            debug=False,
            auto_bold_nss=False,
            t1w_inversion=True,
            name="syn_preprocessing",
        )

        syn_sdc_wf = init_syn_sdc_wf(
            omp_nthreads=threads,
            debug=False,
            name="syn_sdc",
            sloppy=False,
            atlas_threshold=3,
        )

        # read dwi json files to get metadata for each volume
        in_meta = []

        # the b0 is already an average, so we just need the first metadata
        with open(json_file) as f:
            in_meta.append(json.load(f))

        syn_preprocessing_wf.inputs.inputnode.in_epis = Path(in_epis).resolve()
        syn_preprocessing_wf.inputs.inputnode.in_anat = Path(in_anat).resolve()
        syn_preprocessing_wf.inputs.inputnode.std2anat_xfm = Path(
            std2anat_xfm
        ).resolve()
        syn_preprocessing_wf.inputs.inputnode.mask_anat = Path(
            mask_anat
        ).resolve()

        epi_shape = nib.load(in_epis).get_fdata().shape

        if len(epi_shape) == 3:
            nvols = 1
        elif len(epi_shape) == 4:
            nvols = epi_shape[3]

        # create t_masks variable,
        # list of booleans (all true in this case, since all are b0s),
        # length as number of epis
        t_masks = [True for i in range(nvols)]
        syn_preprocessing_wf.inputs.inputnode.t_masks = t_masks
        syn_preprocessing_wf.inputs.inputnode.in_meta = in_meta

        workflow.connect(
            syn_preprocessing_wf,
            "outputnode.epi_ref",
            syn_sdc_wf,
            "inputnode.epi_ref",
        )
        workflow.connect(
            syn_preprocessing_wf,
            "outputnode.anat_ref",
            syn_sdc_wf,
            "inputnode.anat_ref",
        )
        workflow.connect(
            syn_preprocessing_wf,
            "outputnode.anat_mask",
            syn_sdc_wf,
            "inputnode.anat_mask",
        )
        workflow.connect(
            syn_preprocessing_wf,
            "outputnode.sd_prior",
            syn_sdc_wf,
            "inputnode.sd_prior",
        )

        syn_sdc_wf.inputs.inputnode.epi_mask = Path(epi_mask).resolve()

        workflow.base_dir = str(Path(tmpdirname).resolve() / "nipype_wf")

        datasink = pe.Node(nio.DataSink(), name="sinker")
        datasink.inputs.base_directory = str(
            Path(tmpdirname).resolve() / "datasink"
        )

        # mapping from outputnode of syn_sdc_wf to the smk rule output names
        out_file_mapping = {
            "fmap_ref": "unwarped",
            "fmap_mask": "unwarped_mask",
            "fmap": "fmap",
            "out_warp": "xfm",
        }

        for in_name, out_name in out_file_mapping.items():
            workflow.connect(
                syn_sdc_wf, f"outputnode.{in_name}", datasink, f"{out_name}"
            )

        workflow.run()

        # copy file from each datasink folder (sometimes 'plumb' would be in
        # filename, sometimes not, so using glob)
        for name in smk_out.keys():
            copyfile(glob(f"{tmpdirname}/datasink/{name}/*")[0], smk_out[name])


if __name__ == "__main__":
    sdcflows_syn(
        in_epis=snakemake.input.in_epis,  # noqa: F821
        in_anat=snakemake.input.in_anat,  # noqa: F821
        mask_anat=snakemake.input.mask_anat,  # noqa: F821
        std2anat_xfm=snakemake.input.std2anat_xfm,  # noqa: F821
        epi_mask=snakemake.input.epi_mask,  # noqa: F821
        json_file=snakemake.input.json_file,  # noqa: F821
        smk_out=snakemake.output,  # noqa: F821
        threads=snakemake.threads,  # noqa: F821
    )
