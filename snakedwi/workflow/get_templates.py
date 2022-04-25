import os

os.environ["TEMPLATEFLOW_HOME"] = "resources"

from templateflow import api as tflow

str(tflow.get("MNI152NLin2009cSym"))
