rule get_templateflow_files:
    params:
        templateflow_home="resources",
    output:
        t1=config["template_t1w"],
        mask=config["template_mask"],
    container:
        config["singularity"]["python"]
    run:
        import os

        from templateflow import api as tflow

        os.environ["TEMPLATEFLOW_HOME"] = params.templateflow_home
        str(tflow.get(wildcards.template))
