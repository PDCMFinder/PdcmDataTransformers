import os
import pandas as pd
from download_utilities import download_file

def get_CNA_file(row):
    model_name = row[0]
    file_uuid = row[4]
    url = f'https://api.gdc.cancer.gov/data/{file_uuid}'
    file_out = f'../cna/HCMI_cna_{model_name}.tsv'
    os.path.exists(os.path.dirname(file_out))
    print(f'Downloading data from {url}')
    download_file(url, file_out)

# transform CNA count to fold ratios
def transform_CNA_total_count_to_fold_ratio(df):
    df['cna_fold_ratio'] = df['cna_total_count'] / df['cna_total_count']
    return df

case_by_file_list = pd.read_csv("./resources/case-by-file-list.tsv", sep='\t')
CNA_files_index = case_by_file_list['3'].str.contains('gene_level_copy_number.tsv')
case_by_file_list[CNA_files_index].apply(get_CNA_file, axis=1)