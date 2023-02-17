#!/bin/bash
set -euo pipefail
export IFS=$'\n'$'\t'' '

print_help () {
    cat <<EOF
USAGE: mrview_capture img [slice]

Outputs arguments that can be passed to mrview to capture a series of images
across the axial, coronal, and sagital slices

Args
    img     the image to be captured. Used to calculate slice dimensions
    slice   Optional argument specifying the range of captures to be made.
            Follows pythonic START:STOP:STEP syntax. Defaults to ::10
EOF
}

params=()
while [[ -n "${1:-}" && ! "$1" == "--" ]]; do
    case "${1:-help}" in
        --help | -h)
            print_help
            exit 0
            ;;
        -* | --*)
            echo "Error: Unsupported flag $1" >&2
            exit 1
            ;;
        * )
            params["${#params[@]}"]="$1"
            ;;
    esac
    shift
done

if [[ "${#params[@]}" -gt 2 ]]; then
    echo "Error: unrecognized args '${params[@]:2}'" >&2
    exit 1
fi

if [[ "${#params[@]}" -eq 0 ]]; then
    echo "Error: must provide path to image" >&2
    exit 1
fi


# parse and validate image
img="${params[0]}"
if [[ ! -e $img ]]; then
    echo "Error: Image does not exist: '$img'" >&2
    exit 1
fi


# parse and validate slice
eval $( echo "${params[1]:-"::"}" | awk -F':' '
    {print "start="$1" stop="$2" step="$3" err="$4}
')
start="${start:-0}"
step="${step:-10}"

if [[ -n $err ]]; then
    echo 'Error: slice must not have more than 3 numbers (START:STOP:SLICE)' >&2
    exit 1
fi

is_number () {
    [[ "$1" =~ ^-?[0-9]+$ ]]
}

for ix in "$start" "$stop" "$step"; do
    is_number "$ix" || (echo "Error: '$ix' in slice '${params[1]}' is not a number" && exit 1)
done

if [[ "$start" -gt "$stop" ]]; then
    echo "Error: Start value must not be greater than stop: '$start > $stop'" >&2
    exit 1
fi


planes=(sagittal coronal axial)

shape=($(mrinfo -size $img))

i_mid=$(( shape[0] / 2 ))
j_mid=$(( shape[1] / 2 ))
k_mid=$(( shape[2] / 2 ))

for plane_i in 0 1 2; do
    printf -- '-plane %s -capture.prefix %s ' "$plane_i" "${planes[$plane_i]}"
    stop_="${stop:-$(( shape[$plane_i] + step ))}"
    i="$start"
    while [[ "$i" -lt "$stop_" ]]; do
        coord=($i_mid $j_mid $k_mid)
        coord[$plane_i]="$i"
        coord_str=$(echo ${coord[@]} | paste -s -d',')
        printf -- '-voxel %s -capture.grab ' "$coord_str"
        i=$(( i + step ))
    done
done
