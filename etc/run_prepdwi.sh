#!/bin/bash

repo=http://github.com/akhanf/prepdwi_smk
branch=hippdwi

set -euf -o pipefail
confirm() {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case "$response" in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            false
            ;;
    esac
}


if [ "$#" -lt 1 ]
then
 echo "Usage: $0 <out_folder>"
 exit 0
fi

out_folder=$1

if [ ! -e $out_folder ]
then
git clone $repo  $out_folder
#pushd $out_folder
#git checkout $branch
#popd
fi

#edit config
vi $out_folder/config/config.yml

#run dry-run
pushd $out_folder
snakemake -np 
popd

confirm "Run job on interactive node? [y/N]"

export SINGULARITY_NV=1
salloc  -D $out_folder --time=3:00:00 --cpus-per-task=16 --ntasks=1 --mem=64000 --gres=gpu:t4:1 --account=$CC_GPU_ALLOC srun bash snakemake -j16 --use-singularity --res gpus=1 mem_mb=64000

#job_cmd='regularSubmit -g 1 -j 16core64gb24h snakemake -j16 --use-singularity --singularity-args='"\-\-nv"' --res gpus=1 mem_mb=64000'
#
#echo "Submitting job to cluster, with: "
#echo $job_cmd
#
##confirm with user prior to submitting
#confirm 
#
##submit
#$job_cmd

