import argparse
import os
import requests
import pandas as pd
from io import StringIO
from multiprocessing.dummy import Pool

def downloadAndSave(api_root,endpoint_base,filename_base):
    endpoint = f'{api_root}{endpoint_base}'
    filename = f'./data/UPDOG/JAX/JAX_metadata-{filename_base}.tsv'
    print(f'Downloading data from {endpoint}')
    response = requests.get(endpoint)
    response.raise_for_status()
    if not response.text:
        raise IOError(f'No data retrieved from {endpoint}. Check endpoints.')
    with open(filename, 'w') as metadata_file:
        print(f'└── Writing file to {filename}')
        metadata_file.write(response.text)

def get_molecular_data(row):
    api_root = row['api_root']
    endpoint = row['endpoint'] 
    dirname = row['dirname'] 
    model_id = row[0]
    model_endpoint = f'{api_root}{endpoint}={model_id}'
    print(f'Downloading data from {model_endpoint}')
    response = requests.get(model_endpoint) 
    response.raise_for_status()
    molechar_folder = f'./data/UPDOG/JAX/{dirname}'
    if not os.path.exists(molechar_folder):
        os.mkdir(molechar_folder)
    out_file = f'{molechar_folder}/JAX_{dirname}_{model_id}.tsv'
    with open(out_file, 'w+') as molechar_file:
        print(f'└── Saving to file {out_file}')
        molechar_file.write(response.text)
    post_proccess(out_file, dirname) 

def post_proccess(out_file, dirname):
    try:
        omicDf = pd.read_csv(out_file, sep='\t')
        if (~ (omicDf['genome_assembly'] == "GRCh38")).any():
            print("ERROR NON GRCH38 VALUE FOUND. THIS IS A QUALITY CONTROL ISSUE")
        if dirname == "mut":
            omicDf.insert(13, 'strand',"", allow_duplicates=True)
        omicDf.drop(['genome_assembly'], axis=1).to_csv(out_file, sep='\t', index=False)
    except Exception as e: 
        print(e)

def get_molecular_metadata_samples_and_data(api_root):
    filename = f'./data/UPDOG/JAX/JAX_molecular_metadata-sample.tsv'
    with open(filename, 'a+') as molecular_sample:
        endpoints_base_dict = {"modelsWithVariationData": ("model_variation", "mut"), "modelsWithCNAData":("cna","cna"),"modelsWithExpData": ("exp", "expression")}
        for endpoint_base, endpointAndDirname in endpoints_base_dict.items():
            print(f'Downloading {endpoint_base} data')
            endpoint_url = f'{api_root}{endpoint_base}'
            print(f'Retrieving model IDs from {endpoint_url}')
            response = requests.get(endpoint_url) 
            response.raise_for_status()
            response_text = response.text
            molecular_sample.write(response_text)
            model_id_df = pd.read_csv(StringIO(response_text), sep='\t', usecols=['model_id']) 
            model_id_df['api_root'] = api_root
            model_id_df['endpoint'] = endpointAndDirname[0]
            model_id_df['dirname'] = endpointAndDirname[1]
            model_id_df.drop_duplicates(inplace=True)
            pool = Pool(10) 
            pool.map(get_molecular_data, [row for _, row in model_id_df.iterrows()])
            #model_id_df.apply(get_molecular_data, axis=1,  args=(api_root, endpointAndDirname))


def pull_jax_data():
    api_root = "http://tumor.informatics.jax.org/PDXInfo/TSVData.do?"
    domains = ['patient', 'model_validation', 'sharing', 'patient_sample', 'pdx_model']
    for domain in domains:
        downloadAndSave(api_root, domain, domain)

def main():
    api_root = "http://tumor.informatics.jax.org/PDXInfo/TSVData.do?"
    arguments_options = argparse.ArgumentParser("Pulls all metadata and omic data from JAX classes")
    arguments_options.add_argument("-o", "--omics", help="Pull only the JAX molecular data", action="store_true")
    arguments_options.add_argument("-m", "--metadata", help="Pull only JAX molecular data", action="store_true")
    args = arguments_options.parse_args()
    if not os.path.exists('./data/UPDOG/JAX'):
        os.makedirs('./data/UPDOG/JAX')
    if not args.omics:
        pull_jax_data()
    if not args.metadata:
        get_molecular_metadata_samples_and_data(api_root)
main()
