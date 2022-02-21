from os import listdir, getcwd, rename, makedirs
from os.path import isfile, join, isdir, exists
import logging

log = logging.getLogger(__name__)
logging.basicConfig(filename='extract_xlsx_file_updog.log', encoding='utf-8', level=logging.INFO, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%d/%m/%Y %I:%M %p')
start_dir = getcwd()
home = "E:\\EBI\\Work\\data\\updog"

def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]

def get_files(path):
    return [f for f in listdir(path) if isfile(join(path, f))]

def extract_xlsx(path):
    file_list = get_files(path)
    new_path = path.replace('updog', 'data-submission')
    if not exists(new_path):
        makedirs(new_path)
    files_to_extract = [f for f in file_list if f.endswith('.xlsx')]
    if len(files_to_extract)>0:
        for f in files_to_extract:
            rename(join(path, f), join(new_path, f))
            log.info("File extracted: " + f)
    else:
        log.info("No .xlsx file found.")

if not exists(home.replace('updog', 'data-submission')):
    makedirs(home.replace('updog', 'data-submission'))

for provider in get_dirs(home): ## get_dirs will get the provider dirs in updog
    log.info("Working on provider: "+provider)
    extract_xlsx(join(home, provider)) ## File rename: pdx_models to pdx_model using the provider path