from sdcflows.workflows.fit.syn import init_syn_preprocessing_wf
from sdcflows.workflows.fit.syn import init_syn_sdc_wf
import nibabel as nib
from pathlib import Path
import json
from nipype import Workflow
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from shutil import copyfile
import tempfile
from glob import glob

with tempfile.TemporaryDirectory() as tmpdirname:

    workflow = Workflow(name="nipype_wf")

    syn_preprocessing_wf = init_syn_preprocessing_wf(
        omp_nthreads=snakemake.threads,
        debug=False,
        auto_bold_nss=False,
        t1w_inversion=True,
        name=f"syn_preprocessing",
    )

    syn_sdc_wf = init_syn_sdc_wf(
        omp_nthreads=snakemake.threads,
        debug=False,
        name=f"syn_sdc",
        sloppy=False,
        atlas_threshold=3,
    )

    # read dwi json files to get metadata for each volume
    in_meta = []

    # the b0 is already an average, so we just need the first metadata
    with open(snakemake.input.json_file) as f:
        in_meta.append(json.load(f))

    syn_preprocessing_wf.inputs.inputnode.in_epis = Path(
        snakemake.input.in_epis
    ).resolve()
    syn_preprocessing_wf.inputs.inputnode.in_anat = Path(
        snakemake.input.in_anat
    ).resolve()
    syn_preprocessing_wf.inputs.inputnode.std2anat_xfm = Path(
        snakemake.input.std2anat_xfm
    ).resolve()
    syn_preprocessing_wf.inputs.inputnode.mask_anat = Path(
        snakemake.input.mask_anat
    ).resolve()

    epi_shape = nib.load(snakemake.input.in_epis).get_fdata().shape

    if len(epi_shape) == 3:
        nvols = 1
    elif len(epi_shape) == 4:
        nvols = epi_shape[3]

    # create t_masks variable, list of booleans (all true in this case, since all are b0s), length as number of epis
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

    syn_sdc_wf.inputs.inputnode.epi_mask = Path(snakemake.input.epi_mask).resolve()

    workflow.base_dir = str(Path(tmpdirname).resolve() / "nipype_wf")

    datasink = pe.Node(nio.DataSink(), name="sinker")
    datasink.inputs.base_directory = str(Path(tmpdirname).resolve() / "datasink")

    # mapping from outputnode of syn_sdc_wf to the snakemake rule output names
    out_file_mapping = {
        "fmap_ref": "unwarped",
        "fmap_mask": "unwarped_mask",
        "fmap": "fmap",
        "out_warp": "xfm",
    }

    for in_name, out_name in out_file_mapping.items():
        workflow.connect(syn_sdc_wf, f"outputnode.{in_name}", datasink, f"{out_name}")

    workflow.run()

    # copy file from each datasink folder (sometimes 'plumb' would be in filename, sometimes not, so using glob)
    for name in snakemake.output.keys():
        copyfile(glob(f"{tmpdirname}/datasink/{name}/*")[0], snakemake.output[name])
