rule get_templateflow_files:
    params:
        templateflow_home="resources",
    output:
        t1=config["template_t1w"],
        mask=config["template_mask"],
    container: config['singularity']['python']
    script:
        "../scripts/get_templates.py"
