#!/usr/bin/env python
import json
from typing import List

import nibabel as nib
import numpy as np
from numpy import NDArray


def compute_iou(
    tissue_prob_vol: NDArray,
    tissue_k_seg: NDArray,
    sim_prior: NDArray,
    thresh: float,
) -> None:
    s1 = tissue_prob_vol > thresh
    s2 = tissue_k_seg > thresh
    sim_prior = np.sum(np.logical_and(s1, s2).flat) / np.sum(  # noqa: F841
        np.logical_or(s1, s2).flat
    )


def map_channels_to_tissue(
    tissue_priors: List[str],
    seg_channels: str,
    tissue_segs: List[str],
    mapping_json: str,
    smk_config: dict,
) -> None:
    # Load tissue probabilities, warped from template
    tissue_prob_vol = dict()
    for label, nii in zip(smk_config["tissue_labels"], tissue_priors):
        print(f"Label: {label}, Nifti: {nii}")
        tissue_prob_vol[label] = nib.load(nii).get_fdata()

    # Load k-class tissue segmentation
    tissue_k_seg = nib.load(seg_channels)
    tissue_k_seg_shape = tissue_k_seg.shape

    sim_prior_k = np.zeros(
        [len(smk_config["tissue_labels"]), tissue_k_seg_shape[3]]
    )

    # For each prior, find best fitting channel
    for idx, label in enumerate(smk_config["tissue_labels"]):
        for k in range(tissue_k_seg_shape[3]):
            print(f"Computing overlap of {label} prior and channel {k}...")
            # Compute intersection over union
            compute_iou(
                tissue_prob_vol=tissue_prob_vol[label],
                tissue_k_seg=tissue_k_seg.slicer[:, :, :, k].get_fdata(),
                thresh=0.5,
                sim_prior=sim_prior_k[idx, k],
            )
    print(f"Overlap table:\n\n{sim_prior_k}")

    # Write nii to file
    label_to_k_dict = dict()
    for idx, label in enumerate(smk_config["tissue_labels"]):
        label_to_k_dict[label] = int(np.argmax(sim_prior_k[idx, :]))
        print(
            f"Writing image at channel {label_to_k_dict[label]} to output "
            f"file: {tissue_segs[idx]}"
        )
        nib.save(
            tissue_k_seg.slicer[:, :, :, label_to_k_dict[label]],
            tissue_segs[idx],
        )
    with open(mapping_json, "w") as out_fname:
        json.dump(label_to_k_dict, out_fname, indent=4)


if __name__ == "__main__":
    map_channels_to_tissue(
        tissue_priors=snakemake.input.tissue_priors,  # noqa: F821
        seg_channels=snakemake.input.seg_channels,  # noqa: F821
        tissue_segs=snakemake.output.tissue_segs,  # noqa: F821
        mapping_json=snakemake.output.mapping_json,  # noqa: F821
        smk_config=snakemake.config,  # noqa: F821
    )
