


rule affine_to_template:
    input: 
        flo = bids(root='results',suffix='T1w.nii.gz',desc='preproc',datatype='anat',**config['subj_wildcards']),
        ref = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_t1w'])).format(**wildcards),
    output: 
        warped_subj = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz',space='{template}',desc='affine'),
        xfm_ras = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='subject',to='{template}',desc='affine',type_='ras'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.xfm_ras}'

rule convert_template_xfm_ras2itk:
    input:
        bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='ras'),
    output:
        bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='itk'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'c3d_affine_tool {input}  -oitk {output}'

rule warp_brainmask_from_template_affine:
    input: 
        mask = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_mask'])).format(**wildcards),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='subject',to='{template}',desc='affine',type_='itk'),
    output:
        mask = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='mask.nii.gz',from_='{template}',reg='affine',desc='brain'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell: 'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} '
            ' -t [{input.xfm},1] ' #use inverse xfm (going from template to subject)

rule warp_tissue_probseg_from_template_affine:
    input: 
        probseg = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_tissue_probseg'])).format(**wildcards),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='itk'),
    output:
        probseg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='probseg.nii.gz',label='{tissue}',from_='{template}',reg='{desc}'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.probseg} -o {output.probseg} -r {input.ref} '
            ' -t [{input.xfm},1]' #use inverse xfm (going from template to subject)


rule n4biasfield:
    input: 
        t1 = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
        mask = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='mask.nii.gz',from_='{template}'.format(template=config['template']),reg='affine',desc='brain'),
    output:
        t1 = bids(root='work',datatype='anat',**config['subj_wildcards'],desc='n4', suffix='T1w.nii.gz'),
    threads: 8
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}'


rule mask_template_t1w:
    input:
        t1 = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_t1w'])).format(**wildcards),
        mask = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_mask'])).format(**wildcards),
    output:
        t1 = bids(root='work',datatype='anat', prefix='tpl-{template}/tpl-{template}',desc='masked',suffix='T1w.nii.gz')
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'fslmaths {input.t1} -mas {input.mask} {output}'


rule mask_subject_t1w:
    input:
        t1 = bids(root='work',datatype='anat',**config['subj_wildcards'],desc='n4', suffix='T1w.nii.gz'),
        mask = bids(root='work/seg_t1_brain_tissue',**config['subj_wildcards'],suffix='mask.nii.gz',from_='atropos3seg',desc='brain')
    output:
        t1 = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz',from_='atropos3seg',desc='masked'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'fslmaths {input.t1} -mas {input.mask} {output}'


rule ants_syn_affine_init:
    input: 
        flo = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz',from_='atropos3seg',desc='masked'),
        ref = bids(root='work',datatype='anat', prefix='tpl-{template}/tpl-{template}',desc='masked',suffix='T1w.nii.gz'),
        init_xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='subject',to='{template}',desc='affine',type_='itk'),
    params:
        out_prefix = bids(root='work',datatype='anat',suffix='',from_='subject',to='{template}',**config['subj_wildcards']),
        base_opts = '--write-composite-transform -d {dim} --float 1 '.format(dim=config['ants']['dim']),
        intensity_opts = config['ants']['intensity_opts'],
        init_transform = lambda wildcards, input: '-r {xfm}'.format(xfm=input.init_xfm),
        linear_multires = '-c [{reg_iterations},1e-6,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                reg_iterations = config['ants']['linear']['reg_iterations'],
                                shrink_factors = config['ants']['linear']['shrink_factors'],
                                smoothing_factors = config['ants']['linear']['smoothing_factors']),
        linear_metric = lambda wildcards, input: '-m MI[{template},{target},1,32,Regular,0.25]'.format( template=input.ref,target=input.flo),
        deform_model = '-t {deform_model}'.format(deform_model = config['ants']['deform']['transform_model']),
        deform_multires = '-c [{reg_iterations},1e-9,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                reg_iterations = config['ants']['deform']['reg_iterations'],
                                shrink_factors = config['ants']['deform']['shrink_factors'],
                                smoothing_factors = config['ants']['deform']['smoothing_factors']),
        deform_metric = lambda wildcards, input: '-m {metric}[{template},{target},1,4]'.format(
                                metric=config['ants']['deform']['sim_metric'],
                                template=input.ref, target=input.flo)
    output:
        out_composite = bids(root='work',datatype='anat',suffix='Composite.h5',from_='subject',to='{template}',**config['subj_wildcards']),
        out_inv_composite = bids(root='work',datatype='anat',suffix='InverseComposite.h5',from_='subject',to='{template}',**config['subj_wildcards']),
        warped_flo = bids(root='work',datatype='anat',suffix='T1w.nii.gz',space='{template}',desc='SyN',**config['subj_wildcards']),
    threads: 8
    resources:
        mem_mb = 16000, # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time = 60 # 1 hrs
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsRegistration {params.base_opts} {params.intensity_opts} '
        '{params.init_transform} ' #initial xfm  -- rely on this for affine
    #    '-t Rigid[0.1] {params.linear_metric} {params.linear_multires} ' # rigid registration
    #    '-t Affine[0.1] {params.linear_metric} {params.linear_multires} ' # affine registration
        '{params.deform_model} {params.deform_metric} {params.deform_multires} '  # deformable registration
        '-o [{params.out_prefix},{output.warped_flo}]'
       



rule warp_dseg_from_template:
    input: 
        dseg = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_atlas_dseg'])).format(**wildcards),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
        inv_composite = bids(root='work',datatype='anat',suffix='InverseComposite.h5',from_='subject',to='{template}',**config['subj_wildcards']),
    output:
        dseg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}',reg='SyN'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.dseg} -o {output.dseg} -r {input.ref} '
            ' -t {input.inv_composite} ' #use inverse xfm (going from template to subject)


rule warp_tissue_probseg_from_template:
    input: 
        probseg = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_tissue_probseg'])).format(**wildcards),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
        inv_composite = bids(root='work',datatype='anat',suffix='InverseComposite.h5',from_='subject',to='{template}',**config['subj_wildcards']),
    output:
        probseg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='probseg.nii.gz',label='{tissue}',from_='{template}',reg='SyN'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.probseg} -o {output.probseg} -r {input.ref} '
            ' -t {input.inv_composite} ' #use inverse xfm (going from template to subject)

rule warp_brainmask_from_template:
    input: 
        mask = lambda wildcards: workflow.source_path(os.path.join('..','..',config['template_mask'])).format(**wildcards),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz'),
        inv_composite = bids(root='work',datatype='anat',suffix='InverseComposite.h5',from_='subject',to='{template}',**config['subj_wildcards']),
    output:
        mask = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='mask.nii.gz',from_='{template}',reg='SyN',desc='brain'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} '
            ' -t {input.inv_composite} ' #use inverse xfm (going from template to subject)

rule dilate_brainmask:
    input:
        mask = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='mask.nii.gz',from_='{template}',reg='{desc}',desc='brain'),
    params:
        dil_opt =  ' '.join([ '-dilD' for i in range(config['n_init_mask_dilate'])])
    output:
        mask = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='mask.nii.gz',from_='{template}',reg='{desc}',desc='braindilated'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'fslmaths {input} {params.dil_opt} {output}'


#dilate labels N times to provide more of a fudge factor when assigning GM labels
rule dilate_atlas_labels:
    input:
        dseg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}'),
    params:
        dil_opt =  ' '.join([ '-dilD' for i in range(config['n_atlas_dilate'])])
    output:
        dseg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}',desc='dilated'),
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'fslmaths {input} {params.dil_opt} {output}'

rule resample_mask_to_dwi:
    input:
        mask = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='mask.nii.gz',from_='{template}',reg='SyN',desc='brain'),
        ref = bids(root='work',desc='topup',datatype='dwi',method='jac',**config['subj_wildcards'],suffix='b0.nii.gz'),
    params:
        interpolation = 'NearestNeighbor'
    output:
        mask = bids(root='work',**config['subj_wildcards'],desc='brain',suffix='mask.nii.gz',method='template',from_='{template}',reg='SyN',datatype='dwi'),
    container: config['singularity']['ants']
    group: 'subj'
    shell: 
        'antsApplyTransforms -d 3 --input-image-type 0 --input {input.mask} --reference-image {input.ref}  --interpolation {params.interpolation} --output {output.mask} --verbose'
