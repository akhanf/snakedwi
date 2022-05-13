def get_k_tissue_classes(wildcards):
    if wildcards.subject in config["subject_k_tissue_classes"]:
        return config["subject_k_tissue_classes"][wildcards.subject]
    else:
        return config["default_k_tissue_classes"]


# this performs Atropos with k-means as initialization
rule tissue_seg_kmeans_init:
    input:
        t1=bids(
            root="work/reg_t1_to_template",
            **config["subj_wildcards"],
            desc="n4",
            suffix="T1w.nii.gz"
        ),
        mask=bids(
            root="work/reg_t1_to_template",
            **config["subj_wildcards"],
            suffix="mask.nii.gz",
            from_="{template}".format(template=config["template"]),
            reg="affine",
            desc="braindilated"
        ),
    params:
        k=get_k_tissue_classes,
        posterior_fmt="posteriors_%d.nii.gz",
        posterior_glob="posteriors_*.nii.gz",
    output:
        seg=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
            suffix="dseg.nii.gz",
            desc="atroposKseg"
        ),
        posteriors=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
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
        "Atropos -d 3 -a {input.t1} -i KMeans[{params.k}] -x {input.mask} -o [{output.seg},{params.posterior_fmt}] && "
        "fslmerge -t {output.posteriors} {params.posterior_glob} "
        #merge posteriors into a 4d file (intermediate files will be removed b/c shadow)


rule map_channels_to_tissue:
    input:
        tissue_priors=expand(
            bids(
                root="work/reg_t1_to_template",
                **config["subj_wildcards"],
                suffix="probseg.nii.gz",
                label="{tissue}",
                from_="{template}".format(template=config["template"]),
                reg="affine"
            ),
            tissue=config["tissue_labels"],
            allow_missing=True,
        ),
        seg_channels_4d=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
            suffix="probseg.nii.gz",
            desc="atroposKseg"
        ),
    output:
        mapping_json=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
            suffix="mapping.json",
            desc="atropos3seg"
        ),
        tissue_segs=expand(
            bids(
                root="work/seg_t1_brain_tissue",
                **config["subj_wildcards"],
                suffix="probseg.nii.gz",
                label="{tissue}",
                desc="atropos3seg"
            ),
            tissue=config["tissue_labels"],
            allow_missing=True,
        ),
    group:
        "subj"
    container: config['singularity']['python']
    script:
        "../scripts/map_channels_to_tissue.py"


rule tissue_seg_to_4d:
    input:
        tissue_segs=expand(
            bids(
                root="work/seg_t1_brain_tissue",
                **config["subj_wildcards"],
                suffix="probseg.nii.gz",
                label="{tissue}",
                desc="atropos3seg"
            ),
            tissue=config["tissue_labels"],
            allow_missing=True,
        ),
    output:
        tissue_seg=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
            suffix="probseg.nii.gz",
            desc="atropos3seg"
        ),
    group:
        "subj"
    container:
        config["singularity"]["prepdwi"]
    shell:
        "fslmerge -t {output} {input}"


rule brainmask_from_tissue:
    input:
        tissue_seg=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
            suffix="probseg.nii.gz",
            desc="atropos3seg"
        ),
    params:
        threshold=0.5,
    output:
        mask=bids(
            root="work/seg_t1_brain_tissue",
            **config["subj_wildcards"],
            suffix="mask.nii.gz",
            from_="atropos3seg",
            desc="brain"
        ),
    container:
        config["singularity"]["prepdwi"]
    group:
        "subj"
    shell:
        #max over tissue probs, threshold, binarize, fill holes
        "fslmaths {input} -Tmax -thr {params.threshold} -bin -fillh {output}"


# TODO: make lesion mask from the holes in the brainmask (instead of just filling them..) -- could be a nice way to exclude contrast enhanced vessels
