#first, copy to flattened folder structure
rule import:
    input: multiext(config['in_dwi_prefix'],'.nii.gz','.bval','.bvec','.json')
    output:  multiext(bids(root='results',suffix='dwi',**subj_wildcards,**dwi_wildcards),'.nii.gz','.bval','.bvec','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')
 
rule denoise:
    input: bids(root='results',suffix='dwi.nii.gz',**subj_wildcards,**dwi_wildcards) 
    output: bids(root='results',suffix='dwi.nii.gz',desc='denoise',**subj_wildcards,**dwi_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='denoise.log',**subj_wildcards,**dwi_wildcards)
    shell: 'dwidenoise {input} {output} &> {log}' 

rule combine:
    input: expand(bids(root='results',suffix='dwi.nii.gz',desc='denoise',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True)
    output: bids(root='results',suffix='dwi.nii.gz',desc='combined',**subj_wildcards)
    container: config['singularity']['prepdwi']
    shell: 'mrcat {input} {output}'

