rule import_dwi:
    input: multiext(config['in_dwi_prefix'],'.nii.gz','.bval','.bvec','.json')
    output:  multiext(bids(root='results',suffix='dwi',**subj_wildcards,**dwi_wildcards),\
                    '.nii.gz','.bval','.bvec','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')
 
rule denoise:
    input: bids(root='results',suffix='dwi.nii.gz',**subj_wildcards,**dwi_wildcards) 
    output: bids(root='results',suffix='dwi.nii.gz',desc='denoise',**subj_wildcards,**dwi_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='denoise.log',**subj_wildcards,**dwi_wildcards)
    shell: 'dwidenoise {input} {output} 2> {log}' 

rule cp_sidecars_denoise:
    input: multiext(bids(root='results',suffix='dwi',**subj_wildcards,**dwi_wildcards),\
                    '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='denoise',**subj_wildcards,**dwi_wildcards),\
                    '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')


rule run_unring:
    input: bids(root='results',suffix='dwi.nii.gz',desc='denoise',**subj_wildcards,**dwi_wildcards) 
    output: bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='unring.log',**subj_wildcards,**dwi_wildcards)
    shell: 'unring {input} {output} 2> {log}' 

rule cp_sidecars_unring:
    input: multiext(bids(root='results',suffix='dwi',desc='denoise',**subj_wildcards,**dwi_wildcards),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='unring',**subj_wildcards,**dwi_wildcards),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')



def get_dwiextract_cmd(wildcards, input, output):
    if config['flags']['use_avgb0_topup'] == True:
        return f'dwiextract {input.nii} - -fslgrad {input.bvec} {input.bval} -bzero | mrmath - mean {output} -axis 3'
    else:
        return f'dwiextract {input.nii} {output} -fslgrad {input.bvec} {input.bval} -bzero'

rule extract_bzeros:
    input: 
        nii = bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards), 
        bvec = bids(root='results',suffix='dwi.bvec',desc='unring',**subj_wildcards,**dwi_wildcards), 
        bval = bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards,**dwi_wildcards) 
    params:
        #get the avgb0 if use_avgb0_topup is set to true
        cmd = get_dwiextract_cmd 
    output: bids(root='results',suffix='b0.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards) 
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='extract_bzeros.log',**subj_wildcards,**dwi_wildcards)
    group: 'topup'
    shell: '{params.cmd} &> {log}'



#now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii = bids(root='results',suffix='b0.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards),
        json = bids(root='results',suffix='dwi.json',desc='unring',**subj_wildcards,**dwi_wildcards)
    output:
        phenc_txt = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards,**dwi_wildcards),
    group: 'topup'
    script: '../scripts/get_phase_encode_txt.py'
        


rule concat_phase_encode_txt:
    input:
        phenc_txts = expand(bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards,**dwi_wildcards),\
                            **dwi_dict,allow_missing=True)
    output:
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards)
    group: 'topup'
    shell: 'cat {input} > {output}'

rule concat_bzeros:
    input:
        bzero_niis = expand(bids(root='results',suffix='b0.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards),\
                            **dwi_dict,allow_missing=True),
    output:
        bzero_concat = bids(root='results',suffix='b0.nii.gz',desc='unring',**subj_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='concat_bzeros.log',**subj_wildcards)
    group: 'topup'
    shell: 'mrcat {input} {output} 2> {log}'




rule run_topup:
    input:
        bzero_concat = bids(root='results',suffix='b0.nii.gz',desc='unring',**subj_wildcards),
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards)
    params:
        out_prefix = bids(root='results',suffix='topup',**subj_wildcards),
        config = 'b02b0.cnf' #this config sets the multi-res schedule and other params..
    output:
        bzero_corrected = bids(root='results',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
        fieldmap = bids(root='results',suffix='fmap.nii.gz',desc='topup',**subj_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='topup.log',**subj_wildcards)
    group: 'topup'
    shell: 'topup --imain={input.bzero_concat} --datain={input.phenc_concat} --config={params.config}'
           ' --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v --logout={log} '


    
