/*
 * Handle SRA fetching and preprocessing
 */

process GET_ACCESSIONS {
    input:
    file sra_data
    
    output:
    path "acc_list.csv"
    
    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    pd.Series(df.index).to_csv('acc_list.csv', index=False, header=False)
    """ 
}

process GET_AMPLICON_SCHEME {

    input:
    val sra_accession
    file sra_data    

    output:
    val sra_accession
    path('primer_scheme.txt')


    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    scheme = df.loc['${sra_accession}','amplicon_PCR_primer_scheme']
    
    scheme = str(scheme)

    if scheme == 'nan':
        primer_scheme = 'unknown'
    elif 'QIAseq' in scheme or 'v3' in scheme:
        primer_scheme = 'ARTICv3'
    elif 'V5.3' in scheme:
        primer_scheme = 'ARTICv5.3.2'
    elif 'V4.1' or 'v4.1' in scheme:
        primer_scheme = 'ARTICv4.1'
    elif 'SNAP' in scheme:
        primer_scheme = 'snap_primers'
    else:
        primer_scheme = 'unknown'

    with open('primer_scheme.txt', 'w') as f:
            f.write(primer_scheme)
    """
}


process SRA_PREFETCH {

    disk '1TB'
    input:
    val accession
    path primer_scheme

    output:
    tuple val(accession), path(primer_scheme), path("*")

    script:
    """
    #!/bin/sh
    aws s3 sync s3://sra-pub-run-odp/sra/${accession} ${accession} --no-sign-request
    """
}

process FASTERQ_DUMP {
    errorStrategy 'ignore'
    container { params.profile == "docker" ? "ncbi/sra-tools" : "docker://ncbi/sra-tools" }

    input:
    tuple val(accession), path(primer_scheme), path(sra_data)

    output:
    tuple val(accession), path(primer_scheme), path("*.fastq")

    script:
    """
    #!/bin/sh
    fasterq-dump ./${sra_data}/${accession} --progress --threads 8 --split-files    
    """
}