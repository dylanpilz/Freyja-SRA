process AGGREGATE_VARIANTS {
    errorStrategy 'ignore'
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val variants_outputs
    path baseDir

    output:
    path "aggregate_variants.tsv.gz"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    import os
    import subprocess

    paths_list = []
    for file in os.listdir('${baseDir}/outputs/variants'):
        if 'variants' in file:
            paths_list.append(os.path.join('${baseDir}/outputs/variants', file))

    ct_thresh = 20
    j=0
    for var_path in paths_list:
        df = pd.read_csv(var_path,sep='\\t')
        df = df[df['ALT_DP']>=ct_thresh]
        #drop frame shifts
        sname=var_path.split('/')[-1].split('.')[0]
        df = df[df['ALT'].apply(lambda x: True if (('+' not in x and '-' not in x) or ((len(x)-1)%3==0)) else False)]
        if len(df)==0:
            print('No data in ',var_path)
            continue
        df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
        df['sample'] = sname
        df = df[['mutName','sample','ALT_FREQ','ALT_DP']]
        if j==0:
            dfAll = df.copy()
        else:
            dfAll = pd.concat((dfAll,df), axis=0)
        j+=1

    dfAll = dfAll.set_index(['mutName','sample'])
    
    dfAll.to_csv('aggregate_variants.tsv',sep='\\t')

    subprocess.run(["gzip", "aggregate_variants.tsv"])
    
    """
}

process AGGREGATE_DEMIX {

    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    path baseDir

    script:
    """
    python !{baseDir}/scripts/aggregate_demix.py !baseDir
    """
}

process AGGREGATE_COVARIANTS {
    errorStrategy 'ignore'
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val covariants_outputs
    path baseDir

    output:
    path "aggregate_covariants.tsv.gz"

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import os

    agg_df = pd.DataFrame(columns=['Covariants','Sample', 'Count', 'Max_count', 'Freq', 'Coverage_start', 'Coverage_end'])

    for covar_path in os.listdir('${baseDir}/outputs/covariants'):
        df = pd.read_csv('${baseDir}/outputs/covariants/' + covar_path,sep='\\t')
        sample_name = covar_path.split('/')[-1].split('.')[0]
        df['Sample'] = sample_name
        agg_df = pd.concat((agg_df,df), axis=0)

    agg_df = agg_df.set_index(['Covariants','Sample'])
    agg_df.to_csv('aggregate_covariants.tsv',sep='\\t')
    subprocess.run(["gzip", "aggregate_covariants.tsv"])
    """
}