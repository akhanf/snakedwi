from pathlib import Path
from snakebids import bids, generate_inputs, filter_list, get_wildcard_constraints

# from snakeboost import PipEnv


# writes inputs_config.yml and updates config dict
inputs = generate_inputs(
    bids_dir=config["bids_dir"],
    pybids_inputs=config["pybids_inputs"],
    derivatives=config["derivatives"],
    participant_label=config["participant_label"],
    exclude_participant_label=config["exclude_participant_label"],
    use_bids_inputs=True,
)


# this adds constraints to the bids naming
wildcard_constraints:
    **get_wildcard_constraints(config["pybids_inputs"]),


# ---- end snakebids boilerplate ------------------------------------------------


report: "../workflow/report/workflow.rst"


input_wildcards = inputs.input_wildcards
subj_wildcards = inputs.subj_wildcards
input_zip_lists = inputs.input_zip_lists
input_path = inputs.input_path
# setup pipenvs - all my python rules use the script: directive, so will be some work to use snakeboost for this..
# dwi_env = PipEnv(
#            root=Path('work'),
#            requirements=['workflow/pipenvs/snakedwi.txt'],
# )
