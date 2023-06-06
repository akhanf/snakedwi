#!/usr/bin/env python
import json
import sys
from typing import List, Tuple


def check_eddy_s2v(
    eddy_s2v: bool,
    slspec_txt: bool,
    has_slice_timing: bool,
    json_file: str,
) -> None:
    """Check if slice-to-volume correction can be performed"""
    if eddy_s2v and not slspec_txt:
        if not has_slice_timing:
            print(
                "Eddy slice-to-volume enabled, but 'SliceTiming' not found in "
                f"{json_file}.\n"
                "You must do one of the following:\n"
                " 1. add the SliceTiming field to your dwi JSON files\n"
                " 2. use the '--slspec_txt' option to provide a global slspec "
                "file\n"
                " 3. do not set the --use_eddy_s2v CLI option to disable "
                "slice-to-volume correction"
            )
            sys.exit(1)


def check_pe_direction(json_dwi: dict, default_pe_dir: str, json_file: str) -> None:
    """Check if phase encoding direction is defined"""
    if "PhaseEncodingDirection" not in json_dwi:
        if not default_pe_dir == "":
            print("WARNING: setting default PhaseEncodingDirection")
            json_dwi["PhaseEncodingDirection"] = default_pe_dir
        else:
            if "PhaseEncodingAxis" in json_dwi:
                print(
                    "WARNING: assuming PhaseEncodingDirection from " "PhaseEncodingAxis"
                )
                json_dwi["PhaseEncodingDirection"] = json_dwi["PhaseEncodingAxis"]
            else:
                print(
                    f"ERROR: PhaseEncodingDirection not found in {json_file}."
                    "\nYou must add the PhaseEncodingDirection field to your "
                    "dwi JSON files, or use the "
                    "--default_phase_encoding_direction CLI option"
                )
                sys.exit(1)


def check_echo_spacing(json_dwi: dict, default_echo_spacing: str) -> str:
    """Check if effective echo spacing is defined"""
    if "EffectiveEchoSpacing" in json_dwi:
        eff_echo = json_dwi["EffectiveEchoSpacing"]
    elif "EstimatedEffectiveEchoSpacing" in json_dwi:
        eff_echo = json_dwi["EstimatedEffectiveEchoSpacing"]
    else:
        print("EffectiveEchoSpacing not defined in JSON, using default value")
        eff_echo = default_echo_spacing

    return eff_echo


def get_pe_info(json_dwi: dict, pe_axes: List[str], pe_dirs: List[str]) -> None:
    """Extract phase encoding information"""
    pe_axes.append(json_dwi["PhaseEncodingDirection"][0])

    if json_dwi["PhaseEncodingDirection"][-1] == "-":
        pe_dirs.append("-")
    else:
        pe_dirs.append("+")


def get_pe_indices(
    pe_axes: List[str], pe_dirs: List[str], json_files: List[str]
) -> Tuple[List[int], List[str], List[str]]:
    """Extract indices to be used in processing"""
    # If multiple PE axes use LR/RL if available, otherwise AP otherwise
    if len(set(pe_axes)) > 1:
        lr_indices = [idx for idx, pe in enumerate(pe_axes) if pe == "i"]
        ap_indices = [idx for idx, pe in enumerate(pe_axes) if pe == "j"]

        # If opposite PE directions in LR/RL, otherwise use AP
        use_indices = (
            lr_indices
            if len(set([pe_dirs[idx] for idx in lr_indices])) == 2
            else ap_indices
        )
    else:
        use_indices = [idx for idx in range(len(json_files))]

    # Select indices
    pe_axes = [pe_axes[idx] for idx in use_indices]
    pe_dirs = [pe_dirs[idx] for idx in use_indices]

    return use_indices, pe_axes, pe_dirs


def check_subj_dwi_metadata(
    json_files: List[str],
    smk_config: dict,
) -> None:
    bool_map = {True: "yes", False: "no"}
    pe_axes, pe_dirs = [], []
    eddy_s2v = smk_config["use_eddy_s2v"]
    slspec_txt = smk_config["slspec_txt"]

    for json_file in json_files:
        with open(json_file, "r") as fname:
            json_dwi = json.load(fname)

        has_slice_timing = True if "SliceTiming" in json_dwi else False

        # Slice-to-volume
        check_eddy_s2v(
            eddy_s2v=eddy_s2v,
            slspec_txt=slspec_txt,
            has_slice_timing=has_slice_timing,
            json_file=json_file,
        )
        # Phase encoding direction
        check_pe_direction(
            json_dwi=json_dwi,
            default_pe_dir=smk_config["default_phase_encoding_direction"],
            json_file=json_file,
        )
        # Effective echo spacing (not used)
        eff_echo = check_echo_spacing(  # noqa: F841
            json_dwi=json_dwi,
            default_echo_spacing=smk_config["default_effective_echo_spacing"],
        )
        # Save PE information
        get_pe_info(json_dwi=json_dwi, pe_axes=pe_axes, pe_dirs=pe_dirs)

    # Get PE indices
    use_indices, pe_axes, pe_dirs = get_pe_indices(
        pe_axes=pe_axes, pe_dirs=pe_dirs, json_files=json_files
    )

    # Indices
    settings = {
        "indices": use_indices,
        "eddys2v": bool_map[eddy_s2v],
        "PEaxis": pe_axes[0],
    }
    if smk_config["sdc_method"] == "optimal":
        settings["sdc"] = (
            "topup" if len(set(pe_dirs)) >= 2 else smk_config["sdc_method_alternate"]
        )
    else:
        settings["sdc"] = smk_config["sdc_method"]

    return settings
