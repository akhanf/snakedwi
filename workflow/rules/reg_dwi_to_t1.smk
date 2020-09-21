rule import_t1w:
    input:
        t1w = config['in_t1w_preproc'],
    output:
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
    shell: 'cp -v {input.t1w} {output.t1w}'

rule reg_aladin_b0_to_t1:
    input: 
        t1w = bids(root='results',suffix='T1w.nii.gz',desc='preproc',**subj_wildcards),
        avgb0 =  bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards),
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
        avgb0 =  bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards),
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
        res_txt_orig = bids(root='results',suffix='b0.resolution_mm.txt',desc='topup',method='jac',**subj_wildcards),
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
        brainmask = get_mask_for_eddy,
        xfm_itk = bids(root='results',suffix='xfm.txt',from_='dwi',to='T1w',type_='itk',**subj_wildcards),
    params:
        interpolation = 'NearestNeighbor'
    output:
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
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
        brainmask = bids(root='results',suffix='mask.nii.gz',desc='brain',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards),
    params:
        out_basename = lambda wildcards, output: os.path.join(output.out_folder, 'dti')
    output:
        out_folder = directory(bids(root='results',suffix='dtifit',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards)),
        out_fa = os.path.join(directory(bids(root='results',suffix='dtifit',desc='eddy',space='T1w',res=config['resample_dwi']['resample_scheme'],**subj_wildcards)),'dti_FA.nii.gz')
    container: config['singularity']['prepdwi']
    shell:
        'mkdir -p {output.out_folder} && '
        'dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} --mask={input.brainmask} --out={params.out_basename}'


