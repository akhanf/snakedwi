#just grab the first T1w for now:
rule import_t1:
    input: lambda wildcards: expand(config['input_path']['T1w'],zip,**snakebids.filter_list(config['input_zip_lists']['T1w'],wildcards))[0]
    output: bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz')
    group: 'subj'
    shell: 'cp {input} {output}'

rule n4_t1:
    input: 
        t1 = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
    output:
        t1 = bids(root='results',datatype='anat',**config['subj_wildcards'],desc='preproc', suffix='T1w.nii.gz'),
    threads: 8
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'N4BiasFieldCorrection -d 3 -i {input.t1} -o {output}'


rule reg_aladin_b0_to_t1:
    input: 
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',datatype='anat',**config['subj_wildcards']),
        avgb0 =  bids(root='work',suffix='b0.nii.gz',desc='dwiref',datatype='dwi',**config['subj_wildcards']),
    output: 
        warped_avgb0 = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',datatype='dwi',**config['subj_wildcards']),
        xfm_ras = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='ras',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'reg_aladin -rigOnly -flo {input.avgb0} -ref {input.t1w} -res {output.warped_avgb0} -aff {output.xfm_ras}'


rule qc_reg_dwi_t1:
    input:
        ref = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',datatype='dwi',**config['subj_wildcards']),
        flo =  bids(root='results',suffix='T1w.nii.gz',desc='preproc',datatype='anat',**config['subj_wildcards']),
    output:
        png = report(bids(root='qc',suffix='reg.png',**config['subj_wildcards'],from_='dwiref', to='T1w'),
                caption='../report/reg_dwi_t1.rst',
                category='B0 T1w registration'),
        html = bids(root='qc',subject='{subject}',suffix='reg.html',from_='dwiref', to='T1w'),
    group: 'subj'
    script: '../scripts/vis_regqc.py'


rule convert_xfm_ras2itk:
    input:
        xfm_ras = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='ras',datatype='dwi',**config['subj_wildcards']),
    output:
        xfm_itk = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'c3d_affine_tool {input.xfm_ras}  -oitk {output.xfm_itk}'


rule convert_xfm_ras2fsl:
    input:
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',datatype='anat',**config['subj_wildcards']),
        avgb0 =  bids(root='work',suffix='b0.nii.gz',desc='dwiref',datatype='dwi',**config['subj_wildcards']),
        xfm_ras = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='ras',datatype='dwi',**config['subj_wildcards']),
    output:
        xfm_fsl = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='fsl',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'c3d_affine_tool {input.xfm_ras} -ref {input.t1w} -src {input.avgb0} -ras2fsl -o {output.xfm_fsl}'

#tight crop around b0 after rotating into T1w space
rule create_cropped_ref:
    input:
        warped_avgb0 = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',datatype='dwi',**config['subj_wildcards']),
    output:
        cropped_avgb0 = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'c3d {input} -trim 0vox {output}'


#for later resampling.. 
rule write_nii_resolution_to_txt:
    input: '{prefix}.nii.gz'
    output: '{prefix}.resolution_mm.txt'
    group: 'subj'
    script: '../scripts/write_nii_resolution_to_txt.py'
    

#rules for creating reference image for each resampling scheme -- only the rules that are required will be run..
rule create_cropped_ref_t1_resolution:
    input:
        cropped_avgb0 = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',datatype='dwi',**config['subj_wildcards']),
    output:
        avgb0_crop_resample = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',res='T1w',datatype='dwi',**config['subj_wildcards']),
    group: 'subj'
    shell:
        'cp {input} {output}'

rule create_cropped_ref_dwi_resolution:
    input:
        cropped = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',datatype='dwi',**config['subj_wildcards']),
        res_txt_orig = bids(root='work',suffix='b0.resolution_mm.txt',desc='dwiref',datatype='dwi',**config['subj_wildcards']),
    output:
        resampled = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',res='orig',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'c3d {input.cropped} -resample-mm `cat {input.res_txt_orig}` {output}'


rule create_cropped_ref_custom_resolution:
    input:
        cropped = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',datatype='dwi',**config['subj_wildcards']),
    params:
        resolution = 'x'.join([str(vox) for vox in config['resample_dwi']['resample_mm']]) + 'mm'
    output:
        resampled = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',res='custom',datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'c3d {input} -resample-mm {params.resolution} {output}'

rule resample_dwi_to_t1w:
    input:
        ref = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
        dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',datatype='dwi',**config['subj_wildcards']),
        xfm_itk = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',datatype='dwi',**config['subj_wildcards']),
    params:
        interpolation = 'Linear'
    output:
        dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards'])
    container: config['singularity']['ants']
    resources:
        mem_mb = 32000, #-- this is going to be dependent on size of image.. 
    group: 'subj'
    shell:
        'antsApplyTransforms -d 3 --input-image-type 3 --input {input.dwi} --reference-image {input.ref} --transform {input.xfm_itk} --interpolation {params.interpolation} --output {output.dwi} --verbose '
   
rule resample_brainmask_to_t1w:
    input:
        ref = bids(root='work',suffix='avgb0.nii.gz',space='T1w',desc='dwiref',proc='crop',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
        brainmask = get_mask_for_eddy(),
        xfm_itk = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',datatype='dwi',**config['subj_wildcards']),
    params:
        interpolation = 'NearestNeighbor'
    output:
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
    container: config['singularity']['ants']
    resources:
        mem_mb = 32000, #-- this is going to be dependent on size of image.. 
    group: 'subj'
    shell:
        'antsApplyTransforms -d 3 --input-image-type 0 --input {input.brainmask} --reference-image {input.ref} --transform {input.xfm_itk} --interpolation {params.interpolation} --output {output.brainmask} --verbose'
   


rule rotate_bvecs_to_t1w:
    input:
        bvecs = bids(root='results',suffix='dwi.bvec',desc='eddy',datatype='dwi',**config['subj_wildcards']),
        xfm_fsl = bids(root='work',suffix='xfm.txt',from_='dwi',to='T1w',type_='fsl',datatype='dwi',**config['subj_wildcards']),
        bvals = bids(root='results',suffix='dwi.bval',desc='eddy',datatype='dwi',**config['subj_wildcards'])
    params:
        script = os.path.join(config['snakemake_dir'],'workflow/scripts/rotate_bvecs.sh')
    output:
        bvecs = bids(root='results',suffix='dwi.bvec',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
        bvals = bids(root='results',suffix='dwi.bval',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell: 
        '{params.script} {input.bvecs} {input.xfm_fsl} {output.bvecs} && '
        'cp -v {input.bvals} {output.bvals}'    


#dti fitting on dwi in t1w space
rule dtifit_resampled_t1w:
    input:
        dwi = bids(root='results',suffix='dwi.nii.gz',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
        bvals = bids(root='results',suffix='dwi.bval',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
        bvecs = bids(root='results',suffix='dwi.bvec',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards']),
    params:
        out_basename = lambda wildcards, output: os.path.join(output.out_folder, 'dti')
    output:
        out_folder = directory(bids(root='results',suffix='dtifit',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards'])),
        out_fa = os.path.join(directory(bids(root='results',suffix='dtifit',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],datatype='dwi',**config['subj_wildcards'])),'dti_FA.nii.gz')
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'mkdir -p {output.out_folder} && '
        'dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} --mask={input.brainmask} --out={params.out_basename}'

