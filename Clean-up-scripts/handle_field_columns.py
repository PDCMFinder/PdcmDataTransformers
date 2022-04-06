from os import listdir, getcwd, rename, makedirs
from os.path import isfile, join, isdir, exists
import pandas as pd
import logging

log = logging.getLogger(__name__)
logging.basicConfig(filename='handle_field_column_updog.log', level=logging.INFO, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%d/%m/%Y %I:%M %p')
start_dir = getcwd()
home = "/Users/tushar/test-pdxfinder-data/data/UPDOG"

def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]

def get_files(path):
    return [f for f in listdir(path) if isfile(join(path, f))]

def handle_field_column(path):
    tsv_files = [f for f in get_files(path) if f.endswith('.tsv')]
    if len(tsv_files)>0:
        for f in tsv_files:
            metadata = pd.read_csv(join(path,f), sep='\t', na_values="", low_memory=False)
            if 'Field' in metadata.columns:
                metadata = metadata.loc[metadata.Field.str.startswith('#') != True,].reset_index(drop=True)
                metadata = metadata.drop('Field', axis=1)
                metadata.to_csv(join(path, f), sep='\t', index=False)
                log.info("Field column dropped for: " + f)
            else:
                log.info("No field column in the file: "+f)
    else:
        log.info("No .xlsx file found.")

for provider in get_dirs(home): ## get_dirs will get the provider dirs in updog
    log.info("Working on provider: "+provider)
    handle_field_column(join(home, provider)) ## File rename: pdx_models to pdx_model using the provider path