#!/usr/bin/env python

import argparse
from datetime import timedelta,datetime
import os
import json
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET
import http.client
import urllib.error
import time

us_state_to_abbrev = {
    "Alabama": "AL",
    "Alaska": "AK",
    "Arizona": "AZ",
    "Arkansas": "AR",
    "California": "CA",
    "Colorado": "CO",
    "Connecticut": "CT",
    "Delaware": "DE",
    "Florida": "FL",
    "Georgia": "GA",
    "Hawaii": "HI",
    "Idaho": "ID",
    "Illinois": "IL",
    "Indiana": "IN",
    "Iowa": "IA",
    "Kansas": "KS",
    "Kentucky": "KY",
    "Louisiana": "LA",
    "Maine": "ME",
    "Maryland": "MD",
    "Massachusetts": "MA",
    "Michigan": "MI",
    "Minnesota": "MN",
    "Mississippi": "MS",
    "Missouri": "MO",
    "Montana": "MT",
    "Nebraska": "NE",
    "Nevada": "NV",
    "New Hampshire": "NH",
    "New Jersey": "NJ",
    "New Mexico": "NM",
    "New York": "NY",
    "North Carolina": "NC",
    "North Dakota": "ND",
    "Ohio": "OH",
    "Oklahoma": "OK",
    "Oregon": "OR",
    "Pennsylvania": "PA",
    "Rhode Island": "RI",
    "South Carolina": "SC",
    "South Dakota": "SD",
    "Tennessee": "TN",
    "Texas": "TX",
    "Utah": "UT",
    "Vermont": "VT",
    "Virginia": "VA",
    "Washington": "WA",
    "West Virginia": "WV",
    "Wisconsin": "WI",
    "Wyoming": "WY",
    "District of Columbia": "DC",
    "American Samoa": "AS",
    "Guam": "GU",
    "Northern Mariana Islands": "MP",
    "Puerto Rico": "PR",
    "United States Minor Outlying Islands": "UM",
    "U.S. Virgin Islands": "VI",
}

argparser = argparse.ArgumentParser(description='Fetch most recent SRA metadata')

def isnumber(x):
    try:
        float(x)
        return True
    except:
        return False

def get_metadata():

    date_ranges = [
        ('2023-07-01', '2023-10-01'),
        ('2023-04-01', '2023-07-01'),
        ('2023-01-01', '2023-04-01'),
        ('2022-10-01', '2023-01-01'),
        ('2022-07-01', '2022-10-01'),
        ('2022-04-01', '2022-07-01') 
    ]

    metadata = pd.DataFrame()
    for date_range in date_ranges:
        print(date_range)

        start_dt = datetime.strptime(date_range[0], '%Y-%m-%d').date()
        end_dt = datetime.strptime(date_range[1], '%Y-%m-%d').date()

        delta = timedelta(days=1)
        dates = []

        while start_dt <= end_dt:
            dates.append(start_dt.isoformat() + '[All Fields]')
            start_dt += delta

        date_str = ' OR '.join(dates)

        search_term = f'(wastewater[All Fields] AND\
                         ("Severe acute respiratory syndrome coronavirus 2"[Organism] OR sars-cov-2[All Fields])) AND\
                            ({date_str})'

        Entrez.email = "jolevy@scripps.edu"
        handle = Entrez.esearch(db="sra", idtype='acc', retmax=4000,
                                sort='recently_added',
                                term=search_term) 
        record = Entrez.read(handle)
        handle.close()

        try:
            handle = Entrez.efetch(db="sra", id=record['IdList'], rettype="gb",retmode='text')
        except urllib.error.HTTPError as e:
            # Retry once
            print('HTTPError, retrying')
            time.sleep(10)
            handle = Entrez.efetch(db="sra", id=record['IdList'], rettype="gb",retmode='text')

        try:
            string= handle.read()
        except (http.client.IncompleteRead) as e:
            string = e.partial
        handle.close()

        returned_meta=str(string,'UTF-8')

        with open("data/NCBI_metadata.xml", "w") as f:
            f.write(returned_meta)
        
        root = ET.fromstring(returned_meta)
        allDictVals = {}

        for root0 in root:
            ### pull all sample attributes
            vals = [r.text for r in root0.findall('.//SAMPLE_ATTRIBUTE/')]
            sampExp = [r.text for r in root0.findall('.//EXPERIMENT/IDENTIFIERS/PRIMARY_ID')]
            seq_meta = [r.text for r in root0.findall('.//RUN_SET/RUN/RUN_ATTRIBUTES/RUN_ATTRIBUTE/')]
            sampID =  [r.text for r in root0.findall('.//RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID')]
            if len(sampID)>1:
                print('more than one experiment... add funcs')
            elif len(sampID)==0:
                continue
            else:
                sampID = sampID[0]
            ## write to dictionary form
            dictVals = {vals[i].replace(' ','_'):vals[i+1] for i in range(0,len(vals),2)}
            for i in range(0,len(seq_meta),2):
                dictVals[seq_meta[i].replace(' ','')] = seq_meta[i+1]
            dictVals['experiment_id'] = sampID
            dictVals['SRA_id'] = root0[0].attrib['accession']
            allDictVals[sampID] =dictVals

        metadata = pd.concat([metadata, pd.DataFrame(allDictVals).T], axis=0)

    return metadata

def main():
    metadata = get_metadata()

    # Convert collection date to datetime
    metadata  = metadata[metadata['collection_date'].str.contains('20[0-9]{2}-[0-9]{2}-[0-9]{2}')]
    metadata['collection_date'] = pd.to_datetime(metadata['collection_date'].apply(lambda x: x.split('/')[0] if '/' in x else x))

    metadata = metadata.sort_values(by='collection_date',ascending=False)


    # Filter to USA samples
    metadata = metadata[~metadata['geo_loc_name'].isna()]
    metadata = metadata[metadata['geo_loc_name'].str.contains('USA')]

    # Drop duplicates
    metadata = metadata[~metadata.index.duplicated(keep='first')]
    
    # Since Entrez returns the most recent samples, we need to concatenate the new metadata with the old metadata
    current_metadata = pd.read_csv('data/all_metadata.csv', index_col=0)
    new_metadata = metadata[~metadata.index.isin(current_metadata.index)]
    all_metadata = pd.concat([current_metadata, metadata], axis=0)

    all_metadata['ww_population'] = all_metadata['ww_population'].apply(lambda x: x if isnumber(x) else -1.0)
    all_metadata['ww_surv_target_1_conc'] = all_metadata['ww_surv_target_1_conc'].apply(lambda x: x if isnumber(x) else -1.0)

    all_metadata['ww_population'] = all_metadata['ww_population'].fillna(-1.0)
    all_metadata['ww_surv_target_1_conc'] = all_metadata['ww_surv_target_1_conc'].fillna(-1.0)

    all_metadata['geo_loc_country'] = all_metadata['geo_loc_name'].apply(lambda x: x.split(':')[0].strip())
    all_metadata['geo_loc_region'] = all_metadata['geo_loc_name'].apply(lambda x: x.split(':')[1].strip() if len(x.split(':')) > 1 else '')
    all_metadata['geo_loc_region'] = all_metadata['geo_loc_region'].apply(lambda x: x.split(',')[0].strip() if len(x.split(',')) > 1 else x)

    # Create site_id column
    if 'US Virgin Islands' in all_metadata['geo_loc_region'].unique():
        all_metadata['geo_loc_region'] = all_metadata['geo_loc_region'].replace('US Virgin Islands', 'U.S. Virgin Islands')

    all_metadata['site_id'] = (all_metadata['geo_loc_region'].apply(lambda x : us_state_to_abbrev[x])) + all_metadata['ww_population'].fillna('').astype(str)
    all_metadata['site_id'] = all_metadata['site_id'].apply(lambda x: x.replace('.0', ''))

    all_metadata = all_metadata[~all_metadata['site_id'].isna()]

    finished_samples = [file.split('.')[0] for file in os.listdir('outputs/variants')]

    all_metadata = all_metadata[~all_metadata.index.duplicated(keep='first')]

    samples_to_run = all_metadata.copy()
    samples_to_run = samples_to_run[~samples_to_run.index.isin(finished_samples)]

    samples_to_run = samples_to_run[~samples_to_run['ww_surv_target_1_conc'].isna()]
    samples_to_run = samples_to_run[samples_to_run['ww_surv_target_1_conc'] != -1.0]

    samples_to_run = samples_to_run[~samples_to_run['collection_date'].isna()]
    samples_to_run = samples_to_run[~samples_to_run['ww_population'].isna()]
    samples_to_run['ww_population'] = samples_to_run['ww_population'].astype(str)
    

    samples_to_run['collection_date'] = pd.to_datetime(samples_to_run['collection_date'], format='%Y-%m-%d')
    all_metadata['collection_date'] = pd.to_datetime(all_metadata['collection_date'], format='%Y-%m-%d')

    # Temporary filter for freyja paper
    samples_to_run = samples_to_run[samples_to_run['collection_date'] >='2022-04-01']
    samples_to_run = samples_to_run[samples_to_run['collection_date'] <='2023-10-01']

    print('All samples: ', len(all_metadata))
    print('Newly added samples: ', len(new_metadata))
    print('Processed samples: ', len(finished_samples))
    print('Samples to run: ', len(samples_to_run))

    # Sort both dataframes by collection date
    samples_to_run = samples_to_run.sort_values(by='collection_date', ascending=False)
    all_metadata = all_metadata.sort_values(by='collection_date', ascending=False)

    all_metadata.to_csv('data/all_metadata.csv')
    samples_to_run.to_csv('data/samples_to_run.csv', index=True, header=True)

if __name__ == "__main__":
    main()