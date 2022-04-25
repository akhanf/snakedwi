rule get_templateflow_files:
    params:
        templateflow_home="resources",
    output:
        t1=config["template_t1w"],
        mask=config["template_mask"],
    script:
        "../scripts/get_templates.py"
