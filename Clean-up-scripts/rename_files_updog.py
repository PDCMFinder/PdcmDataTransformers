from os import listdir, chdir, getcwd, rename
from os.path import isfile, join, isdir
import logging

log = logging.getLogger(__name__)
logging.basicConfig(filename='rename_file_updog.log', encoding='utf-8', level=logging.INFO, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%d/%m/%Y %I:%M %p')
start_dir = getcwd()
home = "E:\\EBI\\Work\\data\\updog\\"

def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]

def get_files(path):
    return [f for f in listdir(path) if isfile(join(path, f))]

def rename_pdx_model(path, provider):
    file_list = get_files(path)
    target_file, renamed_file = provider + "_metadata-pdx_models.tsv", provider + "_metadata-pdx_model.tsv" ## File names
    if target_file in file_list:
        rename(join(path, target_file), join(path, renamed_file))
        log.info("File renamed: " + target_file)
    else:
        if renamed_file in file_list:
            log.info("File already renamed.")
        else:
            log.warning(target_file + " file not found")

for provider in get_dirs(home): ## get_dirs will get the provider dirs in updog
    log.info("Working on provider: "+provider)
    rename_pdx_model(join(home, provider), provider) ## File rename: pdx_models to pdx_model using the provider path
