#!/usr/bin/env python3

import argparse
import subprocess
import json
import os
import yaml
import pandas as pd

parser = argparse.ArgumentParser(description='Aggregate demix outputs')
parser.add_argument('basedir', type=str, help='Base directory')
args = parser.parse_args()

def isnumber(x):
    try:
        float(x)
        return True
    except:
        return False

def get_alias_key(lineages_yml=f'{args.basedir}/data/lineages.yml'):
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

def main():

    # Create intermediate tsv
    subprocess.run(["freyja", "aggregate", f"{args.basedir}/outputs/demix/", "--output", "aggregate_demix.tsv"])
    
    # Save to json
    agg_demix = pd.read_csv('aggregate_demix.tsv', sep='\t')

    agg_demix['lin_dict'] = [dict(zip(row['lineages'].split(' '), map(float, row['abundances'].split(' ')))) for _, row in agg_demix.iterrows()]
    agg_demix['lin_dict'] = agg_demix['lin_dict'].apply(merge_collapsed)
    agg_demix['lineages'] = agg_demix['lin_dict'].apply(lambda x: ' '.join(list(x.keys())))
    agg_demix['abundances'] = agg_demix['lin_dict'].apply(lambda x: ' '.join([str(v) for v in list(x.values())]))
    agg_demix.drop('lin_dict', axis=1, inplace=True)

    metadata = pd.read_csv(f'{args.basedir}/data/all_metadata.csv')

    columns = ['accession', 'lineages', 'abundances', 'crumbs', 'collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'site_id', 'coverage']

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
    df['coverage'] = agg_demix['coverage']

    metadata.set_index('Unnamed: 0', inplace=True)

    for col in ['collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'site_id']:
        df[col] = df['accession'].map(metadata[col])


    df = df.rename(columns={'ww_surv_target_1_conc':'viral_load'})

    df.set_index('accession', inplace=True)

    df = df[df['lineages'] != ''] 

    df = df[~df['site_id'].isna()]

    with open(f'{args.basedir}/outputs/aggregate/aggregate_demix.json', 'w') as f:
        for row in df.iterrows():

            json_row = {
                'sra_accession': row[0],
                'lineages': [
                    {'name': lineage, 'abundance': float(abundance), 'crumbs': crumbs} for lineage, abundance, crumbs in zip(row[1]['lineages'].split(' '), row[1]['abundances'].split(' '), row[1]['crumbs'])
                ],
                'collection_date': row[1]['collection_date'],
                'geo_loc_country': row[1]['geo_loc_country'],
                'geo_loc_region': row[1]['geo_loc_region'],
                'ww_population': row[1]['ww_population'],
                'viral_load': row[1]['viral_load'],
                'site_id': row[1]['site_id'],
                'coverage': row[1]['coverage']
            }
            
            json_row = json.dumps(json_row)
            f.write(json_row+'\n')

    subprocess.run(['rm', '-rf', 'aggregate_dir'])

if __name__ == '__main__':
    main()