#!/usr/bin/python3
import pandas as pd
import os
import gzip
import requests
import shutil
import json
from glob import glob
from os.path import exists

def get_gz_maf(row):
    link = row["built_links"]
    name = row["Name"]
    if link != "--":
        response = requests.get(link, stream=True)
        print(f'requesting {link}')
        gz_model_filename = response.headers.get('Content-Disposition').replace("attachment; filename=","")
        gz_model_path = f'./mut/{gz_model_filename}' 
        model_filename  = f'./mut/{name}_{gz_model_filename[:-3]}'
        if response.status_code == 200:
            with open(gz_model_path, 'wb') as f:
                f.write(response.raw.read())
            with gzip.open(gz_model_path, 'rb') as f_in:
                print(f'unzipping {gz_model_filename}')
                with open(model_filename, 'wb') as f_out:
                    shutil.copyfileobj(f_in,f_out)
        else:
            print(f'Error downloading at {link}')

def get_biospecimen_file_uuid(row):
    file_uuid = ""
    model_id = row["Name"]
    query_url = f'https://api.gdc.cancer.gov/v0/all?query=nationwidechildrens.org_biospecimen.{model_id}.xml'
    response = requests.get(query_url)
    if response.status_code == 200:
        try:
            file_info = json.loads(response.text)
            file_uuid = file_info.get("data").get("query").get("hits")[0].get("file_id")
            print(f'{model_id} has biospecimen file: {file_uuid}')
        except:
            print(f'Could not retrieve nationwidechildrens.org_biospecimen.{model_id}.xml')
    return file_uuid

def get_biospecimen_file(row):
    model_id = row["Name"]
    file_uuid = row['biospecimen_uuid']
    if file_uuid: 
        query_url = f'https://api.gdc.cancer.gov/data/{file_uuid}'
        print(f'{query_url}') 
        try: 
            response = requests.get(query_url)
            response.raise_for_status()
            with open(f'./xml/CMP_mut_{model_id}.biospecimen.xml', 'w') as xml_file:
                xml_file.write(response.text)
        except:
            print(f'failed to download biospecimen xml for {model_id} with file_uuid {file_uuid}')
            return False

def generate_molecular_metadata(mafDf):
    platform_components = mafDf[['model_id', 'sample_id', 'genome_assembly', 'instrument_model', 'platform_id']]
    unique_platform_components = platform_components.drop_duplicates() 
    molecular_metadata_sample = unique_platform_components[['model_id', 'sample_id','platform_id']]
    molecular_metadata_platform = unique_platform_components[['platform_id', 'genome_assembly','instrument_model']]
    molecular_metadata_platform['molecular_characterisation_type'] = 'mut'
    molecular_metadata_platform['library_strategy'] = 'wxs'
    return molecular_metadata_sample.to_numpy(), molecular_metadata_platform.to_numpy()

def get_uniq_and_save(listOflists, save_file):
    flattened_list  = [row for listOfRows in listOflists for row in listOfRows]
    pd.DataFrame(flattened_list).drop_duplicates().to_csv(save_file, sep='\t', index=False)

def convert_maf_to_pdcm_format():
    columns_to_use = ["Hugo_Symbol", "NCBI_Build","Chromosome","Start_Position","End_Position","Strand","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Mutation_Status","Score","Sequencer","Tumor_Sample_Barcode","t_depth","t_ref_count","t_alt_count","Allele","Existing_variation"]
    vcf_to_pdcm_map = {"Hugo_Symbol": 'symbol', 'NCBI_Build':'genome_assembly', 'Chromosome':'chromosome','Start_Position':'seq_start_pos','End_Position':'seq_end_pos','Strand':'strand','Reference_Allele':'ref_allele', 'Allele':'alt_allele', 'Existing_variation':'variation_id', 't_alt_count': 'read_depth', 'Sequencer': 'instrument_model', 'Tumor_Sample_Barcode': 'sample_id'}
    molecular_metadata_sample = []
    molecular_metadata_platform = []
    print('beginning conversion from vcf to pdcm')
    for maf_file in glob("./mut/*maf"):
        print(f'parsing file {maf_file}') 
        mafDf  = pd.read_csv(maf_file, sep='\t', comment="#", error_bad_lines=False, usecols=columns_to_use)
        pdcm_table = mafDf.rename(columns=vcf_to_pdcm_map)
        model_id = os.path.basename(maf_file.split("_", 1)[0]) 
        pdcm_table['model_id'] = model_id
        pdcm_table['platform_id'] = 'wxs_'  + '_' + pdcm_table['genome_assembly'] + pdcm_table['instrument_model']
        pdcm_table['allele_frequency'] = ""
        sample_entries, platform_entries = generate_molecular_metadata(pdcm_table)
        molecular_metadata_sample.append(sample_entries)
        molecular_metadata_platform.append(platform_entries)
        final_table = pdcm_table[["sample_id","symbol","read_depth","allele_frequency","chromosome","strand","seq_start_pos","seq_end_pos","ref_allele","alt_allele","variation_id","platform_id"]]
        mut_filename = f'./mut/HCMI_mut_{model_id}.tsv' 
        final_table.to_csv(mut_filename, sep='\t', index=False)
    get_uniq_and_save(molecular_metadata_sample, "HCMI_molecular_metadata_sample.tsv")
    get_uniq_and_save(molecular_metadata_platform, "HCMI_molecular_metadata_platform.tsv")

def get_mut_data(cmp):
    cmp['Link To Masked Somatic MAF'].str.replace("https://portal.gdc.cancer.gov/files/", "https://api.gdc.cancer.gov/data/")
    cmp["built_links"] = cmp['Link To Masked Somatic MAF'].str.replace("https://portal.gdc.cancer.gov/files/", "https://api.gdc.cancer.gov/data/")
    cmp.apply(get_gz_maf, axis=1)

def build_biospecimen_uuid(cmp):
    builtTable = None
    uuid_file_path = "./biospecimen_uuid.tsv"
    if not exists(uuid_file_path):
        cmp['biospecimen_uuid'] = cmp.apply(get_biospecimen_file_uuid, axis=1, result_type="expand")
        builtTable = cmp.copy() 
    else:
        builtTable = pd.read_csv(uuid_file_path, sep='\t')
    return builtTable

def get_required_data():
    cmp = pd.read_csv("model-table.tsv", sep='\t', usecols=['Name','Link To Masked Somatic MAF'])
    get_mut_data(cmp) 
    biospecimen_table = build_biospecimen_uuid(cmp)
    biospecimen_table.apply(get_biospecimen_file,axis=1)

convert_maf_to_pdcm_format()
