from sdcflows.workflows.fit.syn import init_syn_preprocessing_wf
from pathlib import Path
import json
from nipype import Workflow
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio

syn_preprocessing_wf = init_syn_preprocessing_wf(
    omp_nthreads=snakemake.threads,
    debug=True,
    auto_bold_nss=False,
    t1w_inversion=True,
    name=f"syn_preprocessing",
)

# read dwi json files (assuming this is the metadata we need)
in_meta = []
for json_file in snakemake.input.json_files:

    with open(json_file) as f:
        in_meta.append(json.load(f))

syn_preprocessing_wf.inputs.inputnode.in_epis = Path(snakemake.input.in_epis).resolve()
syn_preprocessing_wf.inputs.inputnode.in_anat = Path(snakemake.input.in_anat).resolve()
syn_preprocessing_wf.inputs.inputnode.std2anat_xfm = Path(
    snakemake.input.std2anat_xfm
).resolve()
syn_preprocessing_wf.inputs.inputnode.mask_anat = Path(
    snakemake.input.mask_anat
).resolve()
syn_preprocessing_wf.inputs.inputnode.t_masks = [True, True]
syn_preprocessing_wf.inputs.inputnode.in_meta = in_meta


# datasink = pe.Node(nio.DataSink(), name='sinker')
# datasink.inputs.base_directory = str(Path(snakemake.output.out_dir).resolve())

# syn_preprocessing_wf.connect(workflow.connect(inputnode, 'subject_id', datasink, 'container')
syn_preprocessing_wf.base_dir = Path(snakemake.output.out_dir).resolve()

# workflow.connect(syn_preprocessing_wf,'outputnode.epi_ref',datasink,'prep.epi_ref')

syn_preprocessing_wf.run()
# workflow.run()
# init_syn_preprocessing_wf(in_epis=snakemake.input.in_epis,
#                        in_anat=snakemake.input.in_anat,
