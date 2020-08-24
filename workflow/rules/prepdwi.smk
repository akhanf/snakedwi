
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
        fieldmap = bids(root='results',suffix='fmap.nii.gz',desc='topup',**subj_wildcards),
        topup_fieldcoef = bids(root='results',suffix='topup_fieldcoef.nii.gz',**subj_wildcards),
        topup_movpar = bids(root='results',suffix='topup_movpar.txt',**subj_wildcards),

    container: config['singularity']['prepdwi']
    log: bids(root='logs',suffix='topup.log',**subj_wildcards)
    group: 'topup'
    shell: 'topup --imain={input.bzero_concat} --datain={input.phenc_concat} --config={params.config}'
           ' --out={params.out_prefix} --iout={output.bzero_corrected} --fout={output.fieldmap} -v 2> {log}'

#this is for equal positive and negative blipped data -- need to implement the alternative scenario (e.g. b0-only blipped)
rule apply_topup_pos_neg:
    input:
        dwi_niis = expand(bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards,**dwi_wildcards), **dwi_dict, allow_missing=True),
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards),
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
        dwi_topup = bids(root='results',suffix='dwi.nii.gz',desc='topup',**subj_wildcards)
    container: config['singularity']['prepdwi']
    shadow: 'minimal'
    group: 'topup'
    shell: 'applytopup --datain={input.phenc_concat} --imain={params.imain} --inindex={params.inindex} '
           ' -t {params.topup_prefix} -o {params.out_prefix} && '
           ' fslmaths {params.out_prefix}.nii.gz {output.dwi_topup}'



rule cp_sidecars_topup_pos_neg:
    input: multiext(bids(root='results',suffix='dwi',desc='unring',**subj_wildcards,**dwi_exemplar_dict),\
                '.bvec','.bval','.json')
    output: multiext(bids(root='results',suffix='dwi',desc='topup',**subj_wildcards),\
                '.bvec','.bval','.json')
    run:
        for in_file,out_file in zip(input,output):
            shell('cp -v {in_file} {out_file}')


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

       
rule avg_b0_topup:
    input: bids(root='results',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
    output: bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell: 
        'fslmaths {input} -Tmean {output}'

rule cp_brainmask_avg_b0:
    input:
        mask = bids(root='results',subfolder='mask_avgb0',suffix='mask.nii.gz',desc='bet',frac='0.1',smooth='2mm',**subj_wildcards),
    output:
        mask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='avgb0',**subj_wildcards),
    shell:
        'cp -v {input} {output}' 
        
rule cp_brainmask_multishell:
    input:
        mask = bids(root='results',desc='topup',suffix='dwi.avgshells/',**subj_wildcards) + 
                    'atropos_k-6_initmasking_label-brain_smooth-2mm_mask.nii.gz' 
    output:
        mask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='multishell',**subj_wildcards),
    shell:
        'cp -v {input} {output}' 
 
rule qc_brainmask_multishell:
    input: 
        img = bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
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
        img = bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
        seg = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='avgb0',**subj_wildcards),
    output:
        png = report(bids(root='qc',suffix='mask.png',desc='avgb0',**subj_wildcards),
            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False)),

        html = report(bids(root='qc',suffix='mask.html',desc='bet',frac='0.1',smooth='2mm',**subj_wildcards),
            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False))
    script: '../scripts/vis_qc_dseg.py'
        
rule get_slspec_txt:
    input:
        dwi_jsons = expand(bids(root='results',suffix='dwi.json',desc='unring',**subj_wildcards,**dwi_wildcards),**dwi_dict,allow_missing=True),
    output:
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='unring',**subj_wildcards),
    script: '../scripts/get_slspec_txt.py'
         
 
rule run_eddy:
    input:        
        dwi_concat = bids(root='results',suffix='dwi.nii.gz',desc='unring',**subj_wildcards),
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='unring',**subj_wildcards),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='unring',**subj_wildcards),
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='multishell',**subj_wildcards),
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

rule eddy_quad:
    input:
        phenc_concat = bids(root='results',suffix='b0.phenc.txt',desc='unring',**subj_wildcards),
        eddy_index_txt = bids(root='results',suffix='dwi.eddy_index.txt',desc='unring',**subj_wildcards),
        eddy_slspec_txt = bids(root='results',suffix='dwi.eddy_slspec.txt',desc='unring',**subj_wildcards),
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='multishell',**subj_wildcards),
        bvals = bids(root='results',suffix='dwi.bval',desc='unring',**subj_wildcards),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='unring',**subj_wildcards),
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
       
       
#----- T1 registration

rule import_t1w:
    input:
        t1w = config['in_t1w_preproc'],
    output:
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
    shell: 'cp -v {input.t1w} {output.t1w}'

rule reg_aladin_b0_to_t1:
    input: 
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
        avgb0 =  bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
    output: 
        warped_avgb0 = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',**subj_wildcards),
        xfm_ras = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='ras',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'reg_aladin -rigOnly -flo {input.avgb0} -ref {input.t1w} -res {output.warped_avgb0} -aff {output.xfm_ras}'

rule convert_xfm_ras2itk:
    input:
        xfm_ras = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='ras',**subj_wildcards),
    output:
        xfm_itk = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d_affine_tool {input.xfm_ras}  -oitk {output.xfm_itk}'


rule convert_xfm_ras2fsl:
    input:
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
        avgb0 =  bids(root='results',suffix='avgb0.nii.gz',desc='topup',**subj_wildcards),
        xfm_ras = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='ras',**subj_wildcards),
    output:
        xfm_fsl = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='fsl',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d_affine_tool {input.xfm_ras} -ref {input.t1w} -src {input.avgb0} -ras2fsl -o {output.xfm_fsl}'

#tight crop around b0 after rotating into T1w space
rule create_cropped_ref:
    input:
        warped_avgb0 = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',**subj_wildcards),
    output:
        cropped_avgb0 = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d {input} -trim 0vox {output}'


#for later resampling.. 
rule write_nii_resolution_to_txt:
    input: '{prefix}.nii.gz'
    output: '{prefix}.resolution_mm.txt'
    script: '../scripts/write_nii_resolution_to_txt.py'
    

#rules for creating reference image for each resampling scheme -- only the rules that are required will be run..
rule create_cropped_ref_t1_resolution:
    input:
        cropped_avgb0 = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',**subj_wildcards),
    output:
        avgb0_crop_resample = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',res='T1w',**subj_wildcards),
    shell:
        'cp {input} {output}'

rule create_cropped_ref_dwi_resolution:
    input:
        cropped = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',**subj_wildcards),
        res_txt_orig = bids(root='results',suffix='avgb0.resolution_mm.txt',desc='topup',**subj_wildcards),
    output:
        resampled = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',res='orig',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d {input.cropped} -resample-mm `cat {input.res_txt_orig}` {output}'


rule create_cropped_ref_custom_resolution:
    input:
        cropped = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',**subj_wildcards),
    params:
        resolution = 'x'.join([str(vox) for vox in config['resample_dwi']['resample_mm']]) + 'mm'
    output:
        resampled = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',res='custom',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d {input} -resample-mm {params.resolution} {output}'

rule resample_dwi_to_t1w:
    input:
        ref = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
        dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',**subj_wildcards),
        xfm_itk = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',**subj_wildcards),
    params:
        interpolation = 'Linear'
    output:
        dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards)
    container: config['singularity']['ants']
    resources:
        mem_mb = 32000, #-- this is going to be dependent on size of image.. 
    shell:
        'antsApplyTransforms -d 3 --input-image-type 3 --input {input.dwi} --reference-image {input.ref} --transform {input.xfm_itk} --interpolation {params.interpolation} --output {output.dwi} --verbose '
   
rule resample_brainmask_to_t1w:
    input:
        ref = bids(root='results',suffix='avgb0.nii.gz',space='T1w',desc='topup',proc='crop',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',from_='multishell',**subj_wildcards),
        xfm_itk = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',**subj_wildcards),
    params:
        interpolation = 'NearestNeighbor'
    output:
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='topup',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
    container: config['singularity']['ants']
    resources:
        mem_mb = 32000, #-- this is going to be dependent on size of image.. 
    shell:
        'antsApplyTransforms -d 3 --input-image-type 0 --input {input.brainmask} --reference-image {input.ref} --transform {input.xfm_itk} --interpolation {params.interpolation} --output {output.brainmask} --verbose'
   


rule rotate_bvecs_to_t1w:
    input:
        bvecs = bids(root='results',suffix='dwi.bvec',desc='eddy',**subj_wildcards),
        xfm_fsl = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='fsl',**subj_wildcards),
        bvals = bids(root='results',suffix='dwi.bval',desc='eddy',**subj_wildcards)
    output:
        bvecs = bids(root='results',suffix='dwi.bvec',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
        bvals = bids(root='results',suffix='dwi.bval',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards)
    container: config['singularity']['prepdwi']
    shell: 
        'workflow/scripts/rotate_bvecs.sh {input.bvecs} {input.xfm_fsl} {output.bvecs} && '
        'cp -v {input.bvals} {output.bvals}'    


#dti fitting on dwi in t1w space
rule dtifit_resampled_t1w:
    input:
        dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
        bvals = bids(root='results',suffix='dwi.bval',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='topup',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
    params:
        out_basename = lambda wildcards, output: os.path.join(output.out_folder, 'dti')
    output:
        out_folder = directory(bids(root='results',suffix='dtifit',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards)),
        out_fa = os.path.join(directory(bids(root='results',suffix='dtifit',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards)),'dti_FA.nii.gz')
    container: config['singularity']['prepdwi']
    shell:
        'mkdir -p {output.out_folder} && '
        'dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} --mask={input.brainmask} --out={params.out_basename}'


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



