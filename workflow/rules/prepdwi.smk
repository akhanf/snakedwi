from snakebids import bids

wildcard_constraints:
    shell = "[0-9]+",



rule import_dwi:
    input: 
        nii = [re.sub('.nii.gz',ext,config['input_path']['dwi']) for ext in ['.nii.gz','.bval','.bvec','.json']]
    output:
        nii = multiext(bids(root='results',suffix='dwi', datatype='dwi',**config['input_wildcards']['dwi']),'.nii.gz','.bval','.bvec','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule dwidenoise:
    input: multiext(bids(root='results',suffix='dwi',datatype='dwi',**config['input_wildcards']['dwi']),\
                    '.nii.gz','.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='denoise',datatype='dwi',**config['input_wildcards']['dwi']),\
                    '.nii.gz','.bvec','.bval','.json')
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='denoise.log',**config['input_wildcards']['dwi'])
    shell: 'dwidenoise {input[0]} {output[0]} 2> {log} && ' 
            'cp {input[1]} {output[1]} && '
            'cp {input[2]} {output[2]} && '
            'cp {input[3]} {output[3]}'


def get_degibbs_inputs (wildcards):
    # if input dwi at least 30 dirs, then grab denoised as input
    # else grab without denoising 
    import numpy as np
    in_dwi_bval = re.sub('.nii.gz','.bval',config['input_path']['dwi'].format(**wildcards))
    bvals = np.loadtxt(in_dwi_bval)
    if bvals.size < 30:
        prefix = bids(root='results',suffix='dwi',datatype='dwi',**wildcards)
    else:
        prefix = bids(root='results',suffix='dwi',datatype='dwi',desc='denoise',**wildcards)
    return multiext(prefix,'.nii.gz','.bvec','.bval','.json')
 
rule mrdegibbs:
    input: get_degibbs_inputs
    output: multiext(bids(root='results',suffix='dwi',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi']),\
                    '.nii.gz','.bvec','.bval','.json')
    container: config['singularity']['prepdwi']
#    log: bids(root='logs',suffix='degibbs.log',**config['input_wildcards']['dwi'])
    shell: 'mrdegibbs {input[0]} {output[0]} && '#2> {log} && ' 
            'cp {input[1]} {output[1]} && '
            'cp {input[2]} {output[2]} && '
            'cp {input[3]} {output[3]}'


#now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii = bids(root='results',suffix='b0.nii.gz',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi']),
        json = bids(root='results',suffix='dwi.json',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi'])
    output:
        phenc_txt = bids(root='results',suffix='phenc.txt',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi']),
    group: 'topup'
    script: '../scripts/get_phase_encode_txt.py'
        


rule concat_phase_encode_txt:
    input:
        phenc_txts = lambda wildcards: expand(bids(root='results',suffix='phenc.txt',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi']),\
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output:
        phenc_concat = bids(root='results',suffix='phenc.txt',datatype='dwi',desc='degibbs',**config['subj_wildcards'])
    group: 'topup'
    shell: 'cat {input} > {output}'

rule concat_bzeros:
    input:
        bzero_niis = lambda wildcards: expand(bids(root='results',suffix='b0.nii.gz',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi']),\
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output:
        bzero_concat = bids(root='results',suffix='concatb0.nii.gz',datatype='dwi',desc='degibbs',**config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='concat_bzeros.log',**config['subj_wildcards'])
    group: 'topup'
    shell: 'mrcat {input} {output} 2> {log}'


rule run_topup:
    input:
        bzero_concat = bids(root='results',suffix='concatb0.nii.gz',datatype='dwi',desc='degibbs',**config['subj_wildcards']),
        phenc_concat = bids(root='results',suffix='phenc.txt',datatype='dwi',desc='degibbs',**config['subj_wildcards'])
    params:
        out_prefix = bids(root='results',suffix='topup',datatype='dwi',**config['subj_wildcards']),
        config = 'b02b0.cnf' #this config sets the multi-res schedule and other params..
    output:
        bzero_corrected = bids(root='results',suffix='concatb0.nii.gz',desc='topup',datatype='dwi',**config['subj_wildcards']),
        fieldmap = bids(root='results',suffix='fmap.nii.gz',desc='topup',datatype='dwi',**config['subj_wildcards']),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',datatype='dwi',**config['subj_wildcards']),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='topup.log',**config['subj_wildcards'])
    group: 'topup'
    shell: 'topup --imain={input.bzero_concat} --datain={input.phenc_concat} --config={params.config}'
           ' --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v 2> {log}'


#this is for equal positive and negative blipped data - method=lsr
rule apply_topup_lsr:
    input:
        dwi_niis = lambda wildcards: expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',datatype='dwi',**config['input_wildcards']['dwi']),\
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards)),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',datatype='dwi',**config['subj_wildcards']),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',datatype='dwi',**config['subj_wildcards']),
    params:
        #create comma-seperated list of dwi nii
        imain = lambda wildcards, input: ','.join(input.dwi_niis), 
        # create comma-sep list of indices 1-N
        inindex = lambda wildcards, input: ','.join([str(i) for i in range(1,len(input.dwi_niis)+1)]), 
        topup_prefix = bids(root='results',suffix='topup',datatype='dwi',**config['subj_wildcards']),
        out_prefix = 'dwi_topup',
    output: 
        dwi_topup = bids(root='results',suffix='dwi.nii.gz',desc='topup',method='lsr',datatype='dwi',**config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    shadow: 'minimal'
    group: 'topup'
    shell: 'applytopup --verbose --datain={input.phenc_concat} --imain={params.imain} --inindex={params.inindex} '
           ' -t {params.topup_prefix} -o {params.out_prefix} && '
           ' fslmaths {params.out_prefix}.nii.gz {output.dwi_topup}'



rule apply_topup_jac:
    input:
        nii = bids(root='results',suffix='dwi.nii.gz',desc='degibbs',datatype='dwi',**config['input_wildcards']['dwi']), 
        phenc_scan = bids(root='results',suffix='phenc.txt',datatype='dwi',desc='degibbs',**config['input_wildcards']['dwi']),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',datatype='dwi',**config['subj_wildcards']),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',datatype='dwi',**config['subj_wildcards']),
    params:
        topup_prefix = bids(root='results',suffix='topup',datatype='dwi',**config['subj_wildcards']),
    output: 
        nii = bids(root='results',suffix='dwi.nii.gz',desc='topup',method='jac',datatype='dwi',**config['input_wildcards']['dwi']), 
    container: config['singularity']['prepdwi']
    shadow: 'minimal'
    group: 'topup'
    shell: 
        'line=`cat {input.phenc_scan}` && inindex=`grep -n "$line" {input.phenc_concat} | cut -f1 -d:` && '
        ' applytopup --verbose --datain={input.phenc_concat} --imain={input.nii} --inindex=$inindex ' 
        ' -t {params.topup_prefix} -o dwi_topup --method=jac && mv dwi_topup.nii.gz {output.nii}'


#topup-corrected data is only used for brainmasking.. 
# here, use the jac method by default (later can decide if lsr approach can be used based on headers)
# with jac approach, the jac images need to be concatenated, then avgshell extracted

"""
rule cp_sidecars_topup_lsr:
    #TODO: BEST WAY TO TO EXEMPLAR DWI? 
    input: multiext(bids(root='results',suffix='dwi',desc='degibbs',datatype='dwi',**config['subj_wildcards'],**dwi_exemplar_dict),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='topup',method='lsr',datatype='dwi',**config['subj_wildcards']),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')
"""

rule cp_sidecars_topup_jac:
    input: multiext(bids(root='results',suffix='dwi',desc='degibbs',datatype='dwi',**config['subj_wildcards']),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='topup',method='jac',datatype='dwi',**config['subj_wildcards']),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule concat_dwi_topup_jac:
    input:
        dwi_niis = lambda wildcards: expand(bids(root='results',suffix='dwi.nii.gz',desc='topup',method='jac',datatype='dwi',**config['input_wildcards']['dwi']),\
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output:
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='topup',method='jac',datatype='dwi',**config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    group: 'topup'
    shell: 'mrcat {input} {output}' 


rule get_eddy_index_txt:
    input:
        dwi_niis = lambda wildcards: expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',datatype='dwi',**config['input_wildcards']['dwi']),\
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output:
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
    group: 'topup'
    script: '../scripts/get_eddy_index_txt.py'
 
rule concat_degibbs_dwi:
    input:
        dwi_niis = lambda wildcards: expand(bids(root='results',suffix='dwi.nii.gz',desc='degibbs',datatype='dwi',**config['input_wildcards']['dwi']),\
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output:
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='degibbs',datatype='dwi',**config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='concat_degibbs_dwi.log',**config['subj_wildcards'])
    group: 'topup'
    shell: 'mrcat {input} {output} 2> {log}' 

rule concat_runs_bvec:
    input:
        lambda wildcards: expand(bids(root='results',suffix='dwi.bvec',desc='{{desc}}',datatype='dwi',**config['input_wildcards']['dwi']),
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output: bids(root='results',suffix='dwi.bvec',desc='{desc}',datatype='dwi',**config['subj_wildcards'])
    script: '../scripts/concat_bv.py' 

rule concat_runs_bval:
    input:
        lambda wildcards: expand(bids(root='results',suffix='dwi.bval',desc='{{desc}}',datatype='dwi',**config['input_wildcards']['dwi']),
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output: bids(root='results',suffix='dwi.bval',desc='{desc}',datatype='dwi',**config['subj_wildcards'])
    script: '../scripts/concat_bv.py' 

#combines json files from multiple scans -- for now as a hack just copying first json over..
rule concat_runs_json:
    input:
        lambda wildcards: expand(bids(root='results',suffix='dwi.json',desc='{{desc}}',datatype='dwi',**config['input_wildcards']['dwi']),
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output: bids(root='results',suffix='dwi.json',desc='{desc}',datatype='dwi',**config['subj_wildcards'])
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

#have multiple brainmasking workflows -- this rule picks the method chosen in the config file
def get_mask_for_eddy(wildcards):

    #first get name of method
    if wildcards.subject in config['masking']['custom']:
        method = config['masking']['custom'][wildcards.subject]
    else:
        method = config['masking']['default_method']

    #then get bids name of file 
    return bids(root='results',suffix='mask.nii.gz',desc='brain',method=method,datatype='dwi',**config['subj_wildcards'])


#generate qc snapshot for brain  mask 
rule qc_brainmask_for_eddy:
    input:
        img = bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',datatype='dwi',**config['subj_wildcards']),
        seg = get_mask_for_eddy
    output:
#        png = bids(root='qc',subject='{subject}',suffix='mask.png',desc='brain'),
        png = report(bids(root='qc',suffix='mask.png',desc='brain',**config['subj_wildcards']),
                caption='../report/brainmask_dwi.rst',
                category='Brainmask'),

        html = bids(root='qc',suffix='mask.html',desc='brain',**config['subj_wildcards']),
#        html = report(bids(root='qc',subject='{subject}',suffix='dseg.html',atlas='{atlas}', from_='{template}'),
#                caption='../reports/segqc.rst',
#                category='Segmentation QC',
#                subcategory='{atlas} Atlas from {template}'),
    script: '../scripts/vis_qc_dseg.py'

    
rule get_slspec_txt:
    input:
        dwi_jsons = lambda wildcards: expand(bids(root='results',suffix='dwi.json',desc='degibbs',datatype='dwi',**config['input_wildcards']['dwi']),
                            zip,**snakebids.filter_list(config['input_zip_lists']['dwi'], wildcards))
    output:
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
    script: '../scripts/get_slspec_txt.py'
         
 
rule run_eddy:
    input:        
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        brainmask = get_mask_for_eddy,
        bvals = bids(root='results',suffix='dwi.bval',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='degibbs',datatype='dwi',**config['subj_wildcards'])
    params:
        #set eddy output prefix to 'dwi' inside the output folder
        out_prefix = lambda wildcards, output: os.path.join(output.out_folder,'dwi'),
        topup_prefix = bids(root='results',suffix='topup',datatype='dwi',**config['subj_wildcards']),
        flags = ' '.join([f'--{key}' for (key,value) in config['eddy']['flags'].items() if value == True ] ),
        options = ' '.join([f'--{key}={value}' for (key,value) in config['eddy']['opts'].items() if value is not None ] ),
        container = config['singularity']['fsl']
    output:
        #eddy creates many files, so write them to a eddy subfolder instead
        out_folder = directory(bids(root='results',suffix='eddy',datatype='dwi',**config['subj_wildcards'])),
        dwi = os.path.join(bids(root='results',suffix='eddy',datatype='dwi',**config['subj_wildcards']),'dwi.nii.gz'),
        bvec = os.path.join(bids(root='results',suffix='eddy',datatype='dwi',**config['subj_wildcards']),'dwi.eddy_rotated_bvecs')
#    container: config['singularity']['fsl']
    threads: 1
    resources:
        gpus = 1,
        time = 240, #6 hours (this is a conservative estimate, may be shorter)
        mem_mb = 32000,
    log: bids(root='logs',suffix='run_eddy.log',**config['subj_wildcards'])
    group: 'eddy'
    shell: 'singularity exec --nv -e {params.container} eddy_cuda9.1 --imain={input.dwi_concat} --mask={input.brainmask} '
            ' --acqp={input.phenc_concat} --index={input.eddy_index_txt} '
            ' --bvecs={input.bvecs} --bvals={input.bvals} --topup={params.topup_prefix} '
            ' --slspec={input.eddy_slspec_txt} ' 
            ' --out={params.out_prefix} '
            ' {params.flags} {params.options}  &> {log}'


rule cp_eddy_outputs:
    input:
        #get nii.gz, bvec, and bval from eddy output
        dwi = os.path.join(bids(root='results',suffix='eddy',datatype='dwi',**config['subj_wildcards']),'dwi.nii.gz'),
        bvec = os.path.join(bids(root='results',suffix='eddy',datatype='dwi',**config['subj_wildcards']),'dwi.eddy_rotated_bvecs'),
        bval = bids(root='results',suffix='dwi.bval',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
    output:
        multiext(bids(root='results',suffix='dwi',desc='eddy',datatype='dwi',**config['subj_wildcards']),'.nii.gz','.bvec','.bval')
    group: 'eddy'
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')

rule eddy_quad:
    input:
        phenc_concat = bids(root='results',suffix='phenc.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        brainmask = get_mask_for_eddy,
        bvals = bids(root='results',suffix='dwi.bval',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='degibbs',datatype='dwi',**config['subj_wildcards']),
        fieldmap = bids(root='results',suffix='fmap.nii.gz',desc='topup',datatype='dwi',**config['subj_wildcards']),
        eddy_dir = bids(root='results',suffix='eddy',datatype='dwi',**config['subj_wildcards'])
    params: 
        eddy_prefix = lambda wildcards, input: os.path.join(input.eddy_dir,'dwi'),
    output:
        out_dir = directory(bids(root='results',suffix='eddy.qc',datatype='dwi',**config['subj_wildcards'])),
        eddy_qc_pdf = bids(root='results',suffix='eddy.qc/qc.pdf',datatype='dwi',**config['subj_wildcards'])
    
    container: config['singularity']['prepdwi']
    shell: 
        'rmdir {output.out_dir} && '
        'eddy_quad {params.eddy_prefix} --eddyIdx={input.eddy_index_txt} --eddyParams={input.phenc_concat} '
        ' --mask={input.brainmask} --bvals={input.bvals} --bvecs={input.bvecs} --output-dir={output.out_dir} '
        '--slspec={input.eddy_slspec_txt} --verbose'
        #' --field={input.fieldmap} ' #this seems to break it..

rule split_eddy_qc_report:
    input:
        eddy_qc_pdf = bids(root='results',suffix='eddy.qc/qc.pdf',datatype='dwi',**config['subj_wildcards'])
    output:
        report(directory(bids(root='results',suffix='eddy.qc_pages',datatype='dwi',**config['subj_wildcards'])),patterns=['{pagenum}.png'],caption="../report/eddy_qc.rst", category="eddy_qc",subcategory=bids(**config['subj_wildcards'],include_subject_dir=False,include_session_dir=False))
    shell:
        'mkdir -p {output} && convert {input} {output}/%02d.png'
        

       
    #TODO: gradient correction (optional step -- only if gradient file is provided).. 
#  gradient_unwarp.py
#  reg_jacobian
#  convertwarp -> change this to wb_command -convert-warpfield  to get itk transforms 
#  applywarp -> change this to antsApplyTransforms


