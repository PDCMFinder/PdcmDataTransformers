import sys
from os import listdir
from os.path import isfile, join, isdir
import numpy as np
import pandas as pd

def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]


def get_files(path):
    return [join(path, f) for f in listdir(path) if isfile(join(path, f)) and f.endswith(".tsv")]


def read_metadata_without_fields(path):
    metadata = pd.read_csv(path, sep='\t', na_values="", low_memory=False)
    if 'Field' in metadata.columns:
        metadata = metadata.loc[metadata.Field.str.startswith('#') != True,].reset_index(drop=True)
        metadata = metadata.drop('Field', axis=1)
    return metadata


def read_metadata_with_fields(path):
    metadata = pd.read_csv(path, sep='\t', na_values="", low_memory=False)
    return metadata


def sort_case_insensitive(sort_list):
    return sorted(sort_list, key=str.casefold)


def get_hugo2ncbi():
    hugo2ncbi = pd.read_csv("../Clean-up-scripts/hugo_to_ncbi.txt", sep='\t')
    hugo2ncbi["NCBI Gene ID"] = hugo2ncbi["NCBI Gene ID"].fillna(0).astype(int)
    hugo2ncbi["NCBI Gene ID(supplied by NCBI)"] = hugo2ncbi["NCBI Gene ID(supplied by NCBI)"].fillna(0).astype(int)

    hugo2ncbi["NCBI Gene ID"] = np.where(hugo2ncbi["NCBI Gene ID"] != 0, hugo2ncbi["NCBI Gene ID"],
                                         hugo2ncbi["NCBI Gene ID(supplied by NCBI)"])
    hugo2ncbi['RefSeq IDs'] = np.where(hugo2ncbi["RefSeq IDs"].fillna('') != '', hugo2ncbi["RefSeq IDs"],
                                       hugo2ncbi["RefSeq(supplied by NCBI)"].fillna(''))

    hugo2ncbi = hugo2ncbi[hugo2ncbi["NCBI Gene ID"] != 0]
    return hugo2ncbi


def get_mappings(home):
    mappings_path = home.replace('data/UPDOG/', '')
    mappings_path = join(mappings_path, 'mapping/diagnosis_mappings.json')
    mappings = pd.read_json(mappings_path)[['mappingValues', 'mappedTermLabel']]
    mappings = pd.concat([mappings, pd.json_normalize(mappings['mappingValues'])], axis=1)
    mappings = mappings.drop('mappingValues', axis=1).sort_values(by='DataSource').reset_index(drop=True)
    df = pd.read_json(
        'https://dev.cancermodels.org/api/model_metadata?select=cancer_system,histology').drop_duplicates().reset_index(
        drop=True)
    mappings = mappings.merge(df, left_on='mappedTermLabel', right_on='histology', how='left').fillna('-')
    return mappings
