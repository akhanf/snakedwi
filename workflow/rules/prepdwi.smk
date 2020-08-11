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



rule extract_avg_bzero:
    input: 
        nii = bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards), 
        bvec = bids(root='results',suffix='dwi.bvec',desc='unring',**subj_wildcards,**dwi_wildcards), 
        bval = bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards,**dwi_wildcards) 
    output: bids(root='results',suffix='b0.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards) 
    container: config['singularity']['prepdwi']
    group: 'topup'
    shell: 'dwiextract {input.nii} - -fslgrad {input.bvec} {input.bval} -bzero  | mrmath - mean {output} -axis 3'



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

rule avg_b0_topup:
    input: bids(root='results',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
    output: bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell: 
        'fslmaths {input} -Tmean {output}'
 
rule gen_brainmask_for_eddy:
    input: bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
    output: 
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='topup',**subj_wildcards),
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='gen_brainmask_for_eddy.log',**subj_wildcards)
    shell: 
        'bet {input} {output} -f 0.1 &> {log}  && fslmaths {output} -bin {output}'
    
rule get_slspec_txt:
    input:
        dwi_jsons = expand(bids(root='results',suffix='dwi.json',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='unring',**subj_wildcards),
    script: '../scripts/get_slspec_txt.py'
         
   
##def get_eddy_flags:
#    for flags in config['eddy']['flags']:
#        
#    [f'--{}' for value in 
 
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
        #set eddy output prefix to 'dwi' inside the output folder
        out_prefix = lambda wildcards, output: os.path.join(output.out_folder,'dwi'),
        topup_prefix = bids(root='results',suffix='topup',**subj_wildcards),
        flags = ' '.join([f'--{key}' for (key,value) in config['eddy']['flags'].items() if value == True ] ),
        options = ' '.join([f'--{key}={value}' for (key,value) in config['eddy']['opts'].items() if value is not None ] )
    output:
        #eddy creates many files, so write them to a eddy subfolder instead
        out_folder = directory(bids(root='results',suffix='eddy',**subj_wildcards)),
        dwi = os.path.join(bids(root='results',suffix='eddy',**subj_wildcards),'dwi.nii.gz'),
        bvec = os.path.join(bids(root='results',suffix='eddy',**subj_wildcards),'dwi.eddy_rotated_bvecs')
    container: config['singularity']['fsl']
    threads: 1
    resources:
        gpus = 1,
        time = 240, #6 hours (this is a conservative estimate, may be shorter)
        mem_mb = 32000,
    log: bids(root='logs',suffix='run_eddy.log',**subj_wildcards)
    group: 'eddy'
    shell: 'eddy_cuda9.1 --imain={input.dwi_concat} --mask={input.brainmask} '
            ' --acqp={input.phenc_concat} --index={input.eddy_index_txt} '
            ' --bvecs={input.bvecs} --bvals={input.bvals} --topup={params.topup_prefix} '
            ' --slspec={input.eddy_slspec_txt} ' 
            ' --out={params.out_prefix} '
            ' {params.flags} {params.options}  &> {log}'


rule cp_eddy_outputs:
    input:
        #get nii.gz, bvec, and bval from eddy output
        dwi = os.path.join(bids(root='results',suffix='eddy',**subj_wildcards),'dwi.nii.gz'),
        bvec = os.path.join(bids(root='results',suffix='eddy',**subj_wildcards),'dwi.eddy_rotated_bvecs'),
        bval = bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards),
    output:
        multiext(bids(root='results',suffix='dwi',desc='eddy',**subj_wildcards),'.nii.gz','.bvec','.bval')
    group: 'eddy'
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')
       

#----- T1 registration


#use top-up avgb0 to perform registration (so it can be done while eddy is running)

rule import_t1w:
    input:
        t1w = config['in_t1w_preproc'],
    output:
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
    shell: 'cp -v {input.t1w} {output.t1w}'

#rule reg_aladin_b0_to_t1:
    

"""
rule ants_linear_b0_to_t1:
    input: 
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
        b0 = bids(root='results',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
    params:
        out_prefix = bids(root='results',suffix='_',from_='dwi',to='T1w',**subj_wildcards),
        base_opts = '-d {dim} --float 1 --verbose 1 --random-seed {random_seed}'.format(dim=config['ants']['dim'],random_seed=config['ants']['random_seed']),
        intensity_opts = config['ants']['intensity_opts'],
        init_translation = lambda wildcards, input: '-r [{template},{target},1]'.format(template=input.t1w,target=input.b0),
        linear_multires = '-c [{reg_iterations},1e-6,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                reg_iterations = config['ants']['linear']['reg_iterations'],
                                shrink_factors = config['ants']['linear']['shrink_factors'],
                                smoothing_factors = config['ants']['linear']['smoothing_factors']),
        linear_metric = lambda wildcards, input: '-m MI[{template},{target},1,32,Regular,0.25]'.format( template=input.t1w,target=input.b0),
    output:
        out_affine = bids(root='results',suffix='_0GenericAffine.mat',from_='dwi',to='T1w',**subj_wildcards),
        warped_b0 = bids(root='results',suffix='b0.nii.gz',space='T1w',**subj_wildcards),
        warped_t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',space='dwi',**subj_wildcards),
    log: bids(root='logs',suffix='ants_linear_b0_to_t1.log',**subj_wildcards)
    threads: 16
    resources:
        mem_mb = 16000, # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time = 60 # 1 hrs
    container: config['singularity']['ants']
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsRegistration {params.base_opts} {params.intensity_opts} '
        '{params.init_translation} ' #initial translation
        '-t Rigid[0.1] {params.linear_metric} {params.linear_multires} ' # rigid registration
        '-t Affine[0.1] {params.linear_metric} {params.linear_multires} ' # affine registration
        '-o [{params.out_prefix},{output.warped_b0},{output.warped_t1w}] &> {log}'
"""



