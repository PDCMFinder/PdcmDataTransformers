from os import listdir, getcwd, rename, makedirs
from os.path import isfile, join, isdir, exists
import pandas as pd
import logging, itertools
import glob

log = logging.getLogger(__name__)
logging.basicConfig(filename='remove_empty_col_rows.log', level=logging.INFO, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%d/%m/%Y %I:%M %p')
start_dir = getcwd()
home = "/Users/tushar/pdx/pdxfinder-data/data/UPDOG"

def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]


def get_tsv(path, tsvfiles):
    for file in glob.glob(path+"/*.tsv"):
        print(file)
        if isdir(file):
            tsvfiles = get_tsv(join(path, file), tsvfiles)
        else:
            tsvfiles.append(file)
    return tsvfiles


def get_files(path, files):
    [files.append(join(path, f)) if isfile(join(path, f)) else get_files(join(path,f), files) for f in listdir(path)]
    return files

def drop_empty_column(path):
    active_templates = "/Users/tushar/pdx/pdxfinder-data/template/active_templates"
    tsv_files = [f for f in get_files(path, []) if f.endswith('.tsv')]
    if len(tsv_files)>0:
        for f in tsv_files:
            if f.__contains__("sampleplatform-data.tsv") or f.__contains__("raw.tsv") or f.__contains__("web.tsv"):
                continue
            elif f.__contains__("molecular"):
                template_path = join(active_templates, "molecular_metadata")
                #if f.__contains__("platform_web.tsv"):
                #   template_name = "molecular_metadata-platform_web.tsv"
                if f.__contains__("platform.tsv"):
                    template_name = "molecular_metadata-platform.tsv"
                elif f.__contains__("sample.tsv"):
                    template_name = "molecular_metadata-sample.tsv"
            elif f.__contains__("cna"):
                template_path = join(active_templates, "cna")
                template_name = "cna_template-sheet.tsv"
            elif f.__contains__("cytogenetics"):
                template_path = join(active_templates, "cytogenetics")
                template_name = "cytogenetics_template-Sheet1.tsv"
            elif f.__contains__("drug"):
                template_path = join(active_templates, "drug")
                template_name = "drugdosing_template-Sheet1.tsv"
            elif f.__contains__("expression"):
                template_path = join(active_templates, "expression")
                template_name = "expression_template-sheet.tsv"
            elif f.__contains__("mut"):
                template_path = join(active_templates, "mut")
                template_name = "mutation_template_internal-Sheet1.tsv"
            elif f.__contains__("treatment"):
                template_path = join(active_templates, "treatment")
                template_name = "patienttreatment_template-Sheet1.tsv"
            else:
                template_path = join(active_templates, "metadata")
                provider = f.split("/")[-1].split("_")[0]
                template_name = f.split("/")[-1].replace(provider+"_", "").replace("metadata", "metadata_template").replace(" ", "").replace("pdx_models", "pdx_model")
                if f.__contains__("cell_model"):
                    template_name = template_name.replace("cell_model", "cell_model-cell_model")

            template = pd.read_csv(join(template_path, template_name), sep='\t', na_values="", low_memory=False)
            metadata = pd.read_csv(f, sep='\t', na_values="", low_memory=False)
            #metadata = metadata.loc[metadata.Field.str.startswith('#') != True,].reset_index(drop=True)
            #metadata = metadata.drop('Field', axis=1)
            #metadata.to_csv(join(path, f), sep='\t', index=False)
            if len(metadata.columns) != len(template.columns):
                log.info("Columns replaced for: " + f +". WITH: "+template_name)
                log.info("\tNumber of columns: " + str(len(metadata.columns)) + "/" + str(len(template.columns)))
    else:
        log.info("No file found.")

for provider in get_dirs(home): ## get_dirs will get the provider dirs in updog
    log.info("Working on provider: "+provider)
    drop_empty_column(join(home, provider)) ## File rename: pdx_models to pdx_model using the provider path
    log.info("\n\n")