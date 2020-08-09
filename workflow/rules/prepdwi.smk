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
           ' --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v 2> {log}'

#note this will not work unless use_avgb0_topup=True 
rule apply_topup:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards), **dwi_dict, allow_missing=True),
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards)
    params:
        #create comma-seperated list of dwi nii
        imain = lambda wildcards, input: ','.join(input.dwi_niis), 
        # create comma-sep list of indices 1-N
        inindex = lambda wildcards, input: ','.join([str(i) for i in range(1,len(input.dwi_niis)+1)]), 
        topup_prefix = bids(root='results',suffix='topup',**subj_wildcards),
    output: 
        dwi_topup = bids(root='results',suffix='dwi.nii.gz',desc='topup',**subj_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='apply_topup.log',**subj_wildcards)
    group: 'topup'
    shell: 'applytopup --datain={input.phenc_concat} --imain={params.imain} --inindex={params.inindex} '
            '-t {params.topup_prefix} -o {output.dwi_topup} 2> {log}'


rule get_eddy_index_txt:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='unring',**subj_wildcards),
    group: 'topup'
    script: '../scripts/get_eddy_index_txt.py'
 
rule concat_dwi_for_eddy:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='concat_dwi_for_eddy.log',**subj_wildcards)
    group: 'topup'
    shell: 'mrcat {input} {output} 2> {log}' 

rule concat_bvec_for_eddy:
    input:
        expand(bids(root='results',suffix='dwi.bvec',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output: bids(root='results',suffix='dwi.bvec',desc='unring',**subj_wildcards)
    script: '../scripts/concat_bv.py' 

rule concat_bval_for_eddy:
    input:
        expand(bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output: bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards)
    script: '../scripts/concat_bv.py' 


rule gen_brainmask_for_eddy:
    input:
        bzero_corrected = bids(root='results',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
    output: 
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='topup',**subj_wildcards),
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='gen_brainmask_for_eddy.log',**subj_wildcards)
    shell: 
        'bet {input} {output} -f 0.1 2> {log}'
    
rule get_slspec_txt:
    input:
        dwi_jsons = expand(bids(root='results',suffix='dwi.json',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='unring',**subj_wildcards),
    script: '../scripts/get_slspec_txt.py'
         
    
rule run_eddy:
    input:        
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards),
        bvecs_concat = bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards),
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='unring',**subj_wildcards),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='unring',**subj_wildcards),
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='topup',**subj_wildcards),
        bvals = bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='unring',**subj_wildcards)
    params:
        topup_prefix = bids(root='results',suffix='topup',**subj_wildcards),
        out_prefix = bids(root='results',suffix='dwi',desc='eddy',**subj_wildcards),
        eddy_opts = '-v --repol --cnr_maps --residuals --data_is_shelled',
        s2v_opts = '--mporder=6 --s2v_niter=5 --s2v_lambda=1 --s2v_interp=trilinear'
    output:
        out_dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',**subj_wildcards),
        rot_bvec = bids(root='results',suffix='dwi.rotated_bvecs',desc='eddy',**subj_wildcards),
        eddy_params = bids(root='results',suffix='dwi.eddy_parameters',desc='eddy',**subj_wildcards),
        eddy_movement_rms = bids(root='results',suffix='dwi.eddy_movement_rms',desc='eddy',**subj_wildcards),
    container: config['singularity']['prepdwi']
    threads: 8
    resources:
        gpus = 1
    log: bids(root='logs',suffix='run_eddy.log',**subj_wildcards)
    shell: 'eddy_cuda9.1 --imain={input.dwi_concat} --mask={input.brainmask} '
            ' --acqp={input.phenc_concat} --index={input.eddy_index_txt} '
            ' --bvecs={input.bvecs} --bvals={input.bvals} --topup={params.topup_prefix} '
            ' --slspec={input.eddy_slspec_txt} {params.s2v_opts} '
            ' --out={params.out_prefix} {params.eddy_opts} &> {log}'
"""        
        
eddy_openmp --imain=$eddy_work/dwi_uncorrected --mask=$brainmask --acqp=$topup_work/pedir.txt --index=$eddy_work/index.txt --bvecs=$eddy_work/dwi_uncorrected.bvec --bvals=$eddy_work/dwi_uncorrected.bval --topup=$topup_work/topup --out=$eddy_dir/dwi -v --repol

"""        
