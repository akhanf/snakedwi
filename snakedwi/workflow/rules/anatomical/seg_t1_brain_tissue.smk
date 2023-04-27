def get_k_tissue_classes(wildcards):
    if wildcards.subject in config["subject_k_tissue_classes"]:
        return config["subject_k_tissue_classes"][wildcards.subject]
    else:
        return config["default_k_tissue_classes"]


# this performs Atropos with k-means as initialization
rule tissue_seg_kmeans_init:
    input:
        t1=bids(
            root=work,
            datatype="reg_t1_to_template",
            desc="n4",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
        mask=bids(
            root=work,
            datatype="reg_t1_to_template",
            from_="{template}".format(template=config["template"]),
            reg="affine",
            desc="braindilated",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
    params:
        k=get_k_tissue_classes,
        posterior_fmt="posteriors_%d.nii.gz",
        posterior_glob="posteriors_*.nii.gz",
    output:
        seg=bids(
            root=work,
            datatype="seg_t1_brain_tissue",
            **subj_wildcards,
            suffix="dseg.nii.gz",
            desc="atroposKseg"
        ),
        posteriors=bids(
            root=work,
            datatype="seg_t1_brain_tissue",
            **subj_wildcards,
            suffix="probseg.nii.gz",
            desc="atroposKseg"
        ),
    shadow:
        "minimal"
    container:
        config["singularity"]["prepdwi"]
    group:
        "subj"
    shell:
        #merge posteriors into a 4d file 
        #(intermediate files will be removed b/c shadow)
        "Atropos -d 3 -a {input.t1} -i KMeans[{params.k}] -x {input.mask} "
        "-o [{output.seg},{params.posterior_fmt}] && "
        "fslmerge -t {output.posteriors} {params.posterior_glob} "


rule map_channels_to_tissue:
    input:
        tissue_priors=expand(
            bids(
                root="work/reg_t1_to_template",
                **subj_wildcards,
                suffix="probseg.nii.gz",
                label="{tissue}",
                from_="{template}".format(template=config["template"]),
                reg="affine"
            ),
            tissue=config["tissue_labels"],
            allow_missing=True,
        ),
        seg_channels_4d=rules.tissue_seg_kmeans_init.output.posteriors,
    output:
        mapping_json=bids(
            root=work,
            datatype="seg_t1_brain_tissue",
            **subj_wildcards,
            suffix="mapping.json",
            desc="atropos3seg"
        ),
        tissue_segs=expand(
            bids(
                root=work,
                datatype="seg_t1_brain_tissue",
                **subj_wildcards,
                suffix="probseg.nii.gz",
                label="{tissue}",
                desc="atropos3seg"
            ),
            tissue=config["tissue_labels"],
            allow_missing=True,
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/anatomical/map_channels_to_tissue.py"


rule tissue_seg_to_4d:
    input:
        tissue_segs=rules.map_channels_to_tissue.output.tissue_segs,
    output:
        tissue_seg=bids(
            root=work,
            datatype="seg_t1_brain_tissue",
            **subj_wildcards,
            suffix="probseg.nii.gz",
            desc="atropos3seg"
        ),
    group:
        "subj"
    container:
        config["singularity"]["fsl"]
    shell:
        "fslmerge -t {output} {input}"


rule brainmask_from_tissue:
    input:
        tissue_seg=rules.tissue_seg_to_4d.output.tissue_seg,
    params:
        threshold=0.5,
    output:
        mask=bids(
            root=work,
            datatype="seg_t1_brain_tissue",
            **subj_wildcards,
            suffix="mask.nii.gz",
            from_="atropos3seg",
            desc="brain"
        ),
    container:
        config["singularity"]["fsl"]
    group:
        "subj"
    shell:
        #max over tissue probs, threshold, binarize, fill holes
        "fslmaths {input} -Tmax -thr {params.threshold} -bin -fillh {output}"


# TODO: make lesion mask from holes in the brainmask (instead filling them.)
# could be a nice way to exclude contrast enhanced vessels
