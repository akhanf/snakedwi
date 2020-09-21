#input: bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards)
#output: bids(root='results',desc='brain',method='bet_from-b0',suffix='mask.nii.gz',**subj_wildcards)


rule import_avg_b0:
    input:  
        bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards),
    output:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
    shell:
        'cp {input} {output}'

#n4
rule n4_avg_b0:
    input:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='topup',**subj_wildcards),
    output:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='n4',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'N4BiasFieldCorrection -i {input} -o {output}'


#rescale intensities, clip off first/last 5% of intensities, then rescale to 0-2000 
rule rescale_avg_b0:
    input:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='n4',**subj_wildcards),
    output:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='rescale',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d -verbose {input} -clip 5% 95% -stretch 0% 99% 0 2000 -o {output}'

rule bet_avg_b0_default_frac:
    input:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='rescale',**subj_wildcards),
    output:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='bet',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'bet {input} {output}'


rule bet_avg_b0_custom_frac:
    input:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='rescale',**subj_wildcards),
    params:
        frac = '0.{frac}'
    output:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='bet',frac='{frac}',**subj_wildcards)
    container: config['singularity']['prepdwi']
    shell:
        'bet {input} {output} -f {params.frac}'



rule binarize_avg_b0_custom_frac:
    input:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='bet',frac='{frac}',**subj_wildcards),
    output:
        bids(root='results',suffix='mask.nii.gz',desc='brain',method='bet_from-b0',frac='{frac}',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d {input} -binarize  -o {output}'

rule binarize_avg_b0:
    input:
        bids(root='results',subfolder='bet_from-b0',suffix='b0.nii.gz',desc='bet',**subj_wildcards),
    output:
        bids(root='results',suffix='mask.nii.gz',desc='brain',method='bet_from-b0',**subj_wildcards),
    container: config['singularity']['prepdwi']
    shell:
        'c3d {input} -binarize  -o {output}'



