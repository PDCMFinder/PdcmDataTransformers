#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import re
import requests
import json

def clean_case_uuid(row):
     sample_uuid = ""
     link_to_clean = row["Link to Model Details"]
     uuid_group  = re.search('[A-Za-z0-9-]{10,}',link_to_clean)
     if uuid_group:
         sample_uuid = uuid_group.group() 
     return sample_uuid


def get_case_json(case_uuid):
     url = f'https://api.gdc.cancer.gov/cases?{case_uuid}&fields=files.file_id,files.file_name,files.access,files.experimental_strategy'
     print(f'Retrieving {url}')
     response = requests.get(url)
     if response.status_code == 200:
         return response.json()
     else:
         print(f'failed to retrieve {url}')
         return None


def parse_response_json(response_json, case_uuid):
     hits = response_json.get('data').get('hits')
     file_list = []
     if hits:
         for hit in hits:
             files = hit.get("files")
             for file in files:
                file_list.append([case_uuid, file.get('file_id'), file.get('file_name'), file.get('access'), file.get('experimental_strategy')])
     return file_list


def filter_files(file_list):
    filtered_file_list = []
    for file in file_list:
        if file[3] == 'open':
            filtered_file_list.append(file)
    return filtered_file_list

def get_case_files(row):
     case_uuid = row['case_uuid']
     filtered_file_list = []
     if case_uuid:
        response_json = get_case_json(case_uuid)
        if response_json:
             file_list = parse_response_json(response_json,case_uuid)
             filtered_file_list = filter_files(file_list)
     return filtered_file_list


model_table  = pd.read_csv("model-table.tsv", sep='\t', usecols=["Name", "Link to Model Details"])
model_table['case_uuid'] = model_table.apply(clean_case_uuid, axis=1)
case_by_file_list = model_table.apply(get_case_files, axis=1)
flattened_case_by_list = []
for case_files in case_by_file_list:
    for file in case_files:
        flattened_case_by_list.append(file)
pd.DataFrame(flattened_case_by_list).to_csv("./resource/case-by-file-list.tsv", sep='\t', index=False)
