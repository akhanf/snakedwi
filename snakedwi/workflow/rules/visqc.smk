
rule qc_reg:
    input:
        ref=config["template_t1w"],
        flo=bids(
            root="work",
            subject="{subject}",
            suffix="T1w.nii.gz",
            space="{template}",
            desc="{desc}",
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                suffix="regqc.png",
                from_="subject",
                to="{template}",
                desc="{desc}",
            ),
            caption="../report/regqc.rst",
            category="Registration QC",
            subcategory="{desc} {template}",
        ),
        html=bids(
            root="qc",
            subject="{subject}",
            suffix="regqc.html",
            from_="subject",
            to="{template}",
            desc="{desc}",
        ),
    group:
        "subj"
    script:
        "../scripts/vis_regqc.py"


rule qc_probseg:
    input:
        img=bids(root="work", subject="{subject}", desc="n4", suffix="T1w.nii.gz"),
        seg4d=bids(
            root="work",
            subject="{subject}",
            suffix="probseg.nii.gz",
            desc="atropos3seg",
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                suffix="probseg.png",
                desc="atropos3seg",
            ),
            caption="../report/segqc.rst",
            category="Segmentation QC",
            subcategory="3-class Tissue Segmentation",
        ),
    group:
        "subj"
    script:
        "../scripts/vis_qc_probseg.py"


rule qc_dseg:
    input:
        img=bids(root="work", subject="{subject}", desc="n4", suffix="T1w.nii.gz"),
        seg=bids(
            root="work",
            subject="{subject}",
            suffix="dseg.nii.gz",
            atlas="{atlas}",
            from_="{template}",
            reg="SyN",
        ),
    output:
        png=report(
            bids(
                root="qc",
                subject="{subject}",
                suffix="dseg.png",
                atlas="{atlas}",
                from_="{template}",
            ),
            caption="../report/segqc.rst",
            category="Segmentation QC",
            subcategory="{atlas} Atlas from {template}",
        ),
        html=bids(
            root="qc",
            subject="{subject}",
            suffix="dseg.html",
            atlas="{atlas}",
            from_="{template}",
        ),
    group:
        "subj"
    script:
        "../scripts/vis_qc_dseg.py"
