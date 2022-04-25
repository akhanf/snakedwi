import os

os.environ["TEMPLATEFLOW_HOME"] = snakemake.params.templateflow_home

from templateflow import api as tflow

str(tflow.get(snakemake.wildcards.template))
str(tflow.get(snakemake.wildcards.template))
