process AGGREGATE_VARIANTS {
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
    val demix_outputs
    path baseDir

    script:
    """
    #!/usr/bin/env python3

    import subprocess
    import json
    import os
    import yaml
    import pandas as pd

    def get_alias_key(lineages_yml='${baseDir}/data/lineages.yml'):
        with open(lineages_yml, 'r') as alias_key:
            lineage_key = yaml.load(alias_key, Loader=yaml.Loader)
        alias_key = dict([(lin['name'], lin['parent']) for lin in lineage_key if 'parent' in lin])
        alias_key.update([(lin['name'], lin['alias']) for lin in lineage_key if lin['name'] != lin['alias']])
        alias_key.update([r for lin in lineage_key for r in \
            [(lin['name'], lin['name'].split('.')[0]), (lin['name'].split('.')[0], lin['alias'])] \
            if (lin['name'] != lin['alias']) and len(lin['name'].split('.')) == 2 ])
        for n in range(4):
            alias_key.update([(alias, '.'.join(alias.split('.')[:-1])) for name, alias in alias_key.items() if not alias in alias_key and len(alias.split('.')) > 1])
        alias_key.update({'A.1': 'A', 'B.1': 'B'})
        return alias_key

    def _crumbs(lin, alias_key):
        return [lin] + ( _crumbs(alias_key[lin], alias_key) if lin in alias_key else [])

    def crumbs(lin, alias_key):
        lin = lin.upper()
        return _crumbs(lin, alias_key) if lin in alias_key else crumbs(lin[:-1], alias_key) if len(lin.split('.')) > 1 else []

    def merge_collapsed(lin_dict):
        new_dict = {}
        for k in lin_dict.keys():
            if '-like' in k:
                true_lin = k.split('-')[0]
                if true_lin in lin_dict:
                    new_dict[true_lin] = lin_dict[k] + lin_dict[true_lin]
                else:
                    new_dict[true_lin] = lin_dict[k]
            elif 'Misc' not in k:
                new_dict[k] = lin_dict[k]
        return new_dict

    

    subprocess.run(["mkdir", "aggregate_dir"])
    for file in os.listdir('${baseDir}/outputs/demix'):
       subprocess.run(["cp", '${baseDir}/outputs/demix/' + file, "aggregate_dir"])

    # Save to tsv
    subprocess.run(["freyja", "aggregate", "aggregate_dir/", "--output", "aggregate_demix.tsv"])

    # Save to json
    agg_demix = pd.read_csv('aggregate_demix.tsv', sep='\\t')


    agg_demix['lin_dict'] = [dict(zip(row['lineages'].split(' '), map(float, row['abundances'].split(' ')))) for _, row in agg_demix.iterrows()]
    agg_demix['lin_dict'] = agg_demix['lin_dict'].apply(merge_collapsed)
    agg_demix['lineages'] = agg_demix['lin_dict'].apply(lambda x: ' '.join(list(x.keys())))
    agg_demix['abundances'] = agg_demix['lin_dict'].apply(lambda x: ' '.join([str(v) for v in list(x.values())]))
    agg_demix.drop('lin_dict', axis=1, inplace=True)

    metadata = pd.read_csv('${baseDir}/data/all_metadata.csv')
    metadata['geo_loc_country'] = metadata['geo_loc_name'].apply(lambda x: x.split(':')[0].strip())
    metadata['geo_loc_region'] = metadata['geo_loc_name'].apply(lambda x: x.split(':')[1].strip() if len(x.split(':')) > 1 else '')

    columns = ['accession', 'lineages', 'abundances', 'crumbs', 'collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'site_id']

    df = pd.DataFrame(columns=columns)

    agg_demix['Unnamed: 0'] = agg_demix['Unnamed: 0'].apply(lambda x: x.split('.')[0])

    # Drop samples that are not in the metadata
    agg_demix = agg_demix[agg_demix['Unnamed: 0'].isin(metadata['Unnamed: 0'])]

    # Get parent lineage for all lineages
    alias_key = get_alias_key()
    
    df['accession'] = agg_demix['Unnamed: 0']
    df['lineages'] = agg_demix['lineages']
    df['crumbs'] = agg_demix['lineages'].apply(lambda x: [';' + ';'.join(crumbs(lin, alias_key)[::-1]) + ';' for lin in x.split(' ')])
    df['abundances'] = agg_demix['abundances']

    for col in ['collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population','ww_surv_target_1_conc', 'site_id']:
        df[col] = [metadata[metadata['Unnamed: 0'] == x][col] for x in df['accession']]

    # Remove rows where 'ww_surv_target_1_conc' is 'not provided' or 'missing'
    df = df[df['ww_surv_target_1_conc'].astype(str) != 'not provided']
    df = df[df['ww_surv_target_1_conc'].astype(str) != 'missing']
    df = df[df['ww_population'].astype(str) != 'missing']

    df = df.rename(columns={'ww_surv_target_1_conc':'viral_load'})

    df.set_index('accession', inplace=True)

    df = df[df['lineages'] != ''] 

    with open('${baseDir}/outputs/aggregate/aggregate_demix.json', 'w') as f:
        for row in df.iterrows():
            if str(row[1]['ww_population'].values[0]) == 'missing':
                row[1]['ww_population'] = -1.0
            json_row = {
                'sra_accession': row[0],
                'lineages': [
                    {'name': lineage, 'abundance': float(abundance), 'crumbs': crumbs} for lineage, abundance, crumbs in zip(row[1]['lineages'].split(' '), row[1]['abundances'].split(' '), row[1]['crumbs'])
                ],
                'collection_date': row[1]['collection_date'].values[0],
                'geo_loc_country': row[1]['geo_loc_country'].values[0],
                'geo_loc_region': row[1]['geo_loc_region'].values[0],
                'ww_population': float(str(row[1]['ww_population'].values[0]).replace('<','').replace('>','')),
                'viral_load': row[1]['viral_load'].values[0],
                'site_id': row[1]['site_id'].values[0]
            }

            if str(json_row['viral_load']).lower() == 'nan' or json_row['viral_load'] == 'not provided' or json_row['viral_load'] == 'missing':
                json_row['viral_load'] = -1.0
            
            if json_row['ww_population'] == None or str(json_row['ww_population']).lower() == 'nan' or str(json_row['ww_population']).lower() == 'missing':
                json_row['ww_population'] = -1.0
            
            json_row = json.dumps(json_row)
            f.write(json_row+'\\n')

    subprocess.run(["rm", "aggregate_demix.tsv"])
    """
}

process AGGREGATE_COVARIANTS {
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