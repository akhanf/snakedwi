wildcard_constraints:
    shell = "[0-9]+",
    desc ="[a-zA-Z0-9]+",
    acq = "[a-zA-Z0-9]+",
    dir = "[a-zA-Z0-9]+",
    run = "[a-zA-Z0-9]+",


rule import_dwi:
    input: multiext(config['in_dwi_prefix'],'.nii.gz','.bval','.bvec','.json')
    output:  multiext(bids(root='results',suffix='dwi',**subj_wildcards,**dwi_wildcards),\
                    '.nii.gz','.bval','.bvec','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule dwidenoise:
    input: multiext(bids(root='results',suffix='dwi',**subj_wildcards,**dwi_wildcards),\
                    '.nii.gz','.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='denoise',**subj_wildcards,**dwi_wildcards),\
                    '.nii.gz','.bvec','.bval','.json')
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='denoise.log',**subj_wildcards,**dwi_wildcards)
    shell: 'dwidenoise {input[0]} {output[0]} 2> {log} && ' 
            'cp {input[1]} {output[1]} && '
            'cp {input[2]} {output[2]} && '
            'cp {input[3]} {output[3]}'


def get_degibbs_inputs (wildcards):
    # if input dwi at least 30 dirs, then grab denoised as input
    # else grab without denoising 
    import numpy as np
    in_dwi_prefix = config['in_dwi_prefix'].format(**wildcards)
    bvals = np.loadtxt(f'{in_dwi_prefix}.bval')
    if bvals.size < 30:
        prefix = bids(root='results',suffix='dwi',**wildcards)
    else:
        prefix = bids(root='results',suffix='dwi',desc='denoise',**wildcards)
    return multiext(prefix,'.nii.gz','.bvec','.bval','.json')
 
rule mrdegibbs:
    input: get_degibbs_inputs
    output: multiext(bids(root='results',suffix='dwi',desc='degibbs',**subj_wildcards,**dwi_wildcards),\
                    '.nii.gz','.bvec','.bval','.json')
    container: config['singularity']['prepdwi']
#    log: bids(root='logs',suffix='degibbs.log',**subj_wildcards,**dwi_wildcards)
    shell: 'mrdegibbs {input[0]} {output[0]} && '#2> {log} && ' 
            'cp {input[1]} {output[1]} && '
            'cp {input[2]} {output[2]} && '
            'cp {input[3]} {output[3]}'

#this seems to not find bzeros in low-bval scans.. 
# instead of this, we can use the get_shell_avg rule.. 
#rule extract_avg_bzero:
#    input: 
#        nii = bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards), 
#        bvec = bids(root='results',suffix='dwi.bvec',desc='degibbs',**subj_wildcards,**dwi_wildcards), 
#        bval = bids(root='results',suffix='dwi.bval',desc='degibbs',**subj_wildcards,**dwi_wildcards) 
#    output: bids(root='results',suffix='b0.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards) 
#    container: config['singularity']['prepdwi']
#    group: 'topup'
#    shell: 'dwiextract {input.nii} - -fslgrad {input.bvec} {input.bval} -bzero  | mrmath - mean {output} -axis 3'




#now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii = bids(root='results',suffix='b0.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards),
        json = bids(root='results',suffix='dwi.json',desc='degibbs',**subj_wildcards,**dwi_wildcards)
    output:
        phenc_txt = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards,**dwi_wildcards),
    group: 'topup'
    script: '../scripts/get_phase_encode_txt.py'
        


rule concat_phase_encode_txt:
    input:
        phenc_txts = expand(bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards,**dwi_wildcards),\
                            **dwi_dict,allow_missing=True)
    output:
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards)
    group: 'topup'
    shell: 'cat {input} > {output}'

rule concat_bzeros:
    input:
        bzero_niis = expand(bids(root='results',suffix='b0.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards),\
                            **dwi_dict,allow_missing=True),
    output:
        bzero_concat = bids(root='results',suffix='concatb0.nii.gz',desc='degibbs',**subj_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='concat_bzeros.log',**subj_wildcards)
    group: 'topup'
    shell: 'mrcat {input} {output} 2> {log}'


rule run_topup:
    input:
        bzero_concat = bids(root='results',suffix='concatb0.nii.gz',desc='degibbs',**subj_wildcards),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards)
    params:
        out_prefix = bids(root='results',suffix='topup',**subj_wildcards),
        config = 'b02b0.cnf' #this config sets the multi-res schedule and other params..
    output:
        bzero_corrected = bids(root='results',suffix='concatb0.nii.gz',desc='topup',**subj_wildcards),
        fieldmap = bids(root='results',suffix='fmap.nii.gz',desc='topup',**subj_wildcards),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',**subj_wildcards),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',**subj_wildcards),

    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='topup.log',**subj_wildcards)
    group: 'topup'
    shell: 'topup --imain={input.bzero_concat} --datain={input.phenc_concat} --config={params.config}'
           ' --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v 2> {log}'

#this is for equal positive and negative blipped data - method=lsr
rule apply_topup_lsr:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards), **dwi_dict, allow_missing=True),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',**subj_wildcards),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',**subj_wildcards),
    params:
        #create comma-seperated list of dwi nii
        imain = lambda wildcards, input: ','.join(input.dwi_niis), 
        # create comma-sep list of indices 1-N
        inindex = lambda wildcards, input: ','.join([str(i) for i in range(1,len(input.dwi_niis)+1)]), 
        topup_prefix = bids(root='results',suffix='topup',**subj_wildcards),
        out_prefix = 'dwi_topup',
    output: 
        dwi_topup = bids(root='results',suffix='dwi.nii.gz',desc='topup',method='lsr',**subj_wildcards)
    container: config['singularity']['prepdwi']
    shadow: 'minimal'
    group: 'topup'
    shell: 'applytopup --verbose --datain={input.phenc_concat} --imain={params.imain} --inindex={params.inindex} '
           ' -t {params.topup_prefix} -o {params.out_prefix} && '
           ' fslmaths {params.out_prefix}.nii.gz {output.dwi_topup}'

#create command for applytopup -- apply for each phenc dir
def get_applytopup_jac_cmds (wildcards, input,output):
    cmds = []
    topup_prefix = bids(root='results',suffix='topup',**subj_wildcards).format(**wildcards)
    for inindex,(imain,out_nii) in enumerate(zip(input.dwi_niis,output.dwi_niis)):
        out_prefix = bids(root='results',suffix='topup',**subj_wildcards).format(**wildcards),
        cmds.append( f'applytopup --verbose --datain={input.phenc_concat} --imain={imain} --inindex={inindex+1}' \
               f' -t {topup_prefix} -o dwi_topup --method=jac') 
        cmds.append(f'mv dwi_topup.nii.gz {out_nii}')

    return ' && '.join(cmds)

#for applying topup to each phenc dir individually (method=jac)
rule apply_topup_jac:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards), **dwi_dict, allow_missing=True),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',**subj_wildcards),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',**subj_wildcards),
    params:
        #create comma-seperated list of dwi nii
        cmds = get_applytopup_jac_cmds
    output: 
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='topup',method='jac',**subj_wildcards,**dwi_wildcards), **dwi_dict, allow_missing=True)
    container: config['singularity']['prepdwi']
    shadow: 'minimal'
    group: 'topup'
    shell: '{params.cmds}'



#topup-corrected data is only used for brainmasking.. 
# here, use the jac method by default (later can decide if lsr approach can be used based on headers)
# with jac approach, the jac images need to be concatenated, then avgshell extracted


rule cp_sidecars_topup_lsr:
    input: multiext(bids(root='results',suffix='dwi',desc='degibbs',**subj_wildcards,**dwi_exemplar_dict),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='topup',method='lsr',**subj_wildcards),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule cp_sidecars_topup_jac:
    input: multiext(bids(root='results',suffix='dwi',desc='degibbs',**subj_wildcards),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='topup',method='jac',**subj_wildcards),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule concat_dwi_topup_jac:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='topup',method='jac',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='topup',method='jac',**subj_wildcards)
    container: config['singularity']['prepdwi']
    group: 'topup'
    shell: 'mrcat {input} {output}' 


rule get_eddy_index_txt:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='degibbs',**subj_wildcards),
    group: 'topup'
    script: '../scripts/get_eddy_index_txt.py'
 
rule concat_degibbs_dwi:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards)
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='concat_degibbs_dwi.log',**subj_wildcards)
    group: 'topup'
    shell: 'mrcat {input} {output} 2> {log}' 

rule concat_runs_bvec:
    input:
        expand(bids(root='results',suffix='dwi.bvec',desc='{desc}',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output: bids(root='results',suffix='dwi.bvec',desc='{desc}',**subj_wildcards)
    script: '../scripts/concat_bv.py' 

rule concat_runs_bval:
    input:
        expand(bids(root='results',suffix='dwi.bval',desc='{desc}',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output: bids(root='results',suffix='dwi.bval',desc='{desc}',**subj_wildcards)
    script: '../scripts/concat_bv.py' 

#combines json files from multiple scans -- for now as a hack just copying first json over..
rule concat_runs_json:
    input:
        expand(bids(root='results',suffix='dwi.json',desc='{desc}',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output: bids(root='results',suffix='dwi.json',desc='{desc}',**subj_wildcards)
    shell: 'cp {input[0]} {output}'
#    script: '../scripts/concat_json.py' 

rule get_shells_from_bvals:
    input: '{dwi_prefix}.bval'
    output: '{dwi_prefix}.shells.json'
    script:
        '../scripts/get_shells_from_bvals.py'
 
#writes 4d file
rule get_shell_avgs:
    input: 
        dwi = '{dwi_prefix}.nii.gz',
        shells = '{dwi_prefix}.shells.json'
    output: 
        avgshells = '{dwi_prefix}.avgshells.nii.gz'
    script:
        '../scripts/get_shell_avgs.py'

#this gets a particular shell (can use to get b0)
rule get_shell_avg:
    input:
        dwi = '{dwi_prefix}_dwi.nii.gz',
        shells = '{dwi_prefix}_dwi.shells.json'
    params:
        bval = '{shell}'
    output:
        avgshell = '{dwi_prefix}_b{shell}.nii.gz'
    script:
        '../scripts/get_shell_avg.py'
        
#unused?       
#rule avg_b0_topup:
#    input: bids(root='results',suffix='concatb0.nii.gz',desc='topup',**subj_wildcards),
#    output: bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
#    container: config['singularity']['prepdwi']
#    shell: 
#        'fslmaths {input} -Tmean {output}'


#have multiple brainmasking workflows -- then a rule to import the one we want..
#this is the output:
#mask = bids(root='results',suffix='mask.nii.gz',desc='brain',for_='eddy',**subj_wildcards),

def get_mask_for_eddy(wildcards):

    #first get name of method
    if wildcards.subject in config['masking']['custom']:
        method = config['masking']['custom'][wildcards.subject]
    else:
        method = config['masking']['default_method']

    #then get bids name of file 
    return bids(root='results',suffix='mask.nii.gz',desc='brain',method=method,**subj_wildcards)

    
"""
rule cp_brainmask_avg_b0:
    input:
        mask = bids(root='results',subfolder='mask_avgb0',suffix='mask.nii.gz',desc='bet',frac='0.1',smooth='2mm',**subj_wildcards),
    output:
        mask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='avgb0',**subj_wildcards),
    shell:
        'cp -v {input} {output}' 
        
rule cp_brainmask_multishell:
    input:
        mask = bids(root='results',desc='topup',method='jac',suffix='dwi.avgshells/',**subj_wildcards) + 
                    'atropos_k-6_initmasking_label-brain_smooth-2mm_mask.nii.gz' 
    output:
        mask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='multishell',**subj_wildcards),
    shell:
        'cp -v {input} {output}' 
 
rule qc_brainmask_multishell:
    input: 
        img = bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards),
        seg = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='multishell',**subj_wildcards),
    output:
        png = report(bids(root='qc',suffix='mask.png',desc='multishell',**subj_wildcards),
            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False)),

        html = report(bids(root='qc',suffix='mask.html',desc='multishell',k='6',smooth='2mm',**subj_wildcards),
            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False))
    script: '../scripts/vis_qc_dseg.py'
 
rule qc_brainmask_avg_b0:
    input: 
        img = bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards),
        seg = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='avgb0',**subj_wildcards),
    output:
        png = report(bids(root='qc',suffix='mask.png',desc='avgb0',**subj_wildcards),
            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False)),

        html = report(bids(root='qc',suffix='mask.html',desc='bet',frac='0.1',smooth='2mm',**subj_wildcards),
            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False))
    script: '../scripts/vis_qc_dseg.py'
   
"""
     
rule get_slspec_txt:
    input:
        dwi_jsons = expand(bids(root='results',suffix='dwi.json',desc='degibbs',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='degibbs',**subj_wildcards),
    script: '../scripts/get_slspec_txt.py'
         
 
rule run_eddy:
    input:        
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='degibbs',**subj_wildcards),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='degibbs',**subj_wildcards),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='degibbs',**subj_wildcards),
        brainmask = get_mask_for_eddy,
        bvals = bids(root='results',suffix='dwi.bval',desc='degibbs',**subj_wildcards),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='degibbs',**subj_wildcards)
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
        bval = bids(root='results',suffix='dwi.bval',desc='degibbs',**subj_wildcards),
    output:
        multiext(bids(root='results',suffix='dwi',desc='eddy',**subj_wildcards),'.nii.gz','.bvec','.bval')
    group: 'eddy'
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule eddy_quad:
    input:
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',**subj_wildcards),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='degibbs',**subj_wildcards),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='degibbs',**subj_wildcards),
        brainmask = get_mask_for_eddy,
        bvals = bids(root='results',suffix='dwi.bval',desc='degibbs',**subj_wildcards),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='degibbs',**subj_wildcards),
        fieldmap = bids(root='results',suffix='fmap.nii.gz',desc='topup',**subj_wildcards),
        eddy_dir = bids(root='results',suffix='eddy',**subj_wildcards)
    params: 
        eddy_prefix = lambda wildcards, input: os.path.join(input.eddy_dir,'dwi'),
    output:
        out_dir = directory(bids(root='results',suffix='eddy.qc',**subj_wildcards)),
        eddy_qc_pdf = bids(root='results',suffix='eddy.qc/qc.pdf',**subj_wildcards)
    
    container: config['singularity']['prepdwi']
    shell: 
        'rmdir {output.out_dir} && '
        'eddy_quad {params.eddy_prefix} --eddyIdx={input.eddy_index_txt} --eddyParams={input.phenc_concat} '
        ' --mask={input.brainmask} --bvals={input.bvals} --bvecs={input.bvecs} --output-dir={output.out_dir} '
        '--slspec={input.eddy_slspec_txt} --verbose'
        #' --field={input.fieldmap} ' #this seems to break it..

rule split_eddy_qc_report:
    input:
        eddy_qc_pdf = bids(root='results',suffix='eddy.qc/qc.pdf',**subj_wildcards)
    output:
        report(directory(bids(root='results',suffix='eddy.qc_pages',**subj_wildcards)),patterns=['{pagenum}.png'],caption="../report/eddy_qc.rst", category="eddy_qc",subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False))
    shell:
        'mkdir -p {output} && convert {input} {output}/%02d.png'
        
#rule eddy_quad_report: #need this separate from eddy_quad, since adding to report seems to create folder too
#    input: 
#        out_dir = bids(root='results',suffix='eddy.qc',**subj_wildcards)
#    output: 
#        qc_pdf = report(bids(root='results',suffix='eddy.qc/qc.pdf',**subj_wildcards))
       
    #TODO: gradient correction (optional step -- only if gradient file is provided).. 
#  gradient_unwarp.py
#  reg_jacobian
#  convertwarp -> change this to wb_command -convert-warpfield  to get itk transforms 
#  applywarp -> change this to antsApplyTransforms




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



