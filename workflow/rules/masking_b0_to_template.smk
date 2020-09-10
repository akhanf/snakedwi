
# should create an avgb0 template in MNI2009cSym space..  - can use ants for that.. single iteration to initial template

rule reg_b0_to_template:
    input:  
        flo = bids(root='results',suffix='b0.nii.gz',desc='topup',method='jac',**subj_wildcards),
        
    output:

