import sys
from os import listdir
from os.path import isfile, join
import pandas as pd
import logging

log = logging.getLogger(__name__)
logging.basicConfig(filename='data_merge.log', level=logging.INFO, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%d/%m/%Y %I:%M %p')

def get_files(path, files):
    [files.append(join(path, f)) if isfile(join(path, f)) else get_files(join(path,f), files) for f in listdir(path)]
    return [f for f in files if f.__contains__('.tsv')]


def remove_field_col(path):
    metadata = pd.read_csv(path, sep='\t', na_values="", low_memory=False)
    if 'Field' in metadata.columns:
        metadata = metadata.loc[metadata.Field.str.startswith('#') != True,].reset_index(drop=True).drop('Field', axis=1)
    return metadata


def data_merge():
    source, target = sys.argv[1], sys.argv[2]
    source_data = get_files(source, [])
    target_data = get_files(target, [])
    target_files_list = [f.replace(target, "") for f in target_data]
    source_files_list = [f.replace(source, "") for f in source_data]
    common_files = list(set(target_files_list) & set(source_files_list))
    for f in common_files:
        source_metadata, target_metadata = remove_field_col(join(source, f)), remove_field_col(join(target, f))
        if len(list(set(target_metadata.iloc[:,0]) - set(source_metadata.iloc[:,0]))) > 0:
            log.info('New data in '+f)
    print(sys.argv)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        data_merge()