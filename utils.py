from os import listdir, makedirs, getcwd
from os.path import isfile, join, isdir, exists
import pandas as pd
from tqdm import tqdm
import requests
from random import randint
import numpy as np


geneSymbol_location = pd.read_csv("/Users/tushar/CancerModels/utils/PdcmDataTransformers/resources/genes.tsv", sep='\t')
home = "/Users/tushar/CancerModels/pdxfinder-data/data/UPDOG/"
templates = "/Users/tushar/CancerModels/pdxfinder-data/template/active_templates"

# Define a function to convert log2 ratio to fold change
def log2_to_fold_change(log2_ratio):
    return 2 ** log2_ratio

# Define a function to convert log2 ratio to log10 ratio
def log2_to_log10(log2_ratio):
    return log2_ratio / np.log10(2)

def request(link, flag, req_type):
    success_quotes = {0: 'May the force be with you! ',
                      1: 'Lights will guide you home! :*',
                      2: 'Voila, got it!',
                      3: 'Look what I found!!!'}

    fail_quotes = {0: 'Oops something went wrong!',
                   1: "You've ran into some trouble, check your input",
                   2: "What I've done!!!! - Linkin Park",
                   3: "Guess what I couldn't find."}
    try:
        if req_type == "get":
            response = requests.get(link)
        elif req_type == "post":
            response = requests.post(link)
        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Parse the JSON data (assuming the API returns JSON)
            if flag:
                print(success_quotes[randint(0, len(success_quotes)-1)])
            return response
        else:
            # If the request was not successful, raise an exception
            if flag:
                print(fail_quotes[randint(0, len(fail_quotes)-1)])
            response.raise_for_status()
    except requests.exceptions.RequestException as e:
        # Handle any exceptions that may occur during the request
        print(f"An error occurred: {e}")
        return None



def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]


def get_files(path):
    return [f for f in listdir(path) if isfile(join(path, f))]


def read_metadata_without_fields(path):
    metadata = pd.read_csv(path, sep='\t', na_values="", low_memory=False, encoding='utf-8', encoding_errors='replace')
    if 'Field' in metadata.columns:
        metadata = metadata.loc[metadata.Field.astype(str).str.startswith('#') != True,].reset_index(drop=True)
        metadata = metadata.drop('Field', axis=1)
    return metadata


def read_metadata_with_fields(path):
    metadata = pd.read_csv(path, sep='\t', na_values="", low_memory=False, encoding='utf-8', encoding_errors='replace')
    return metadata


def sort_case_insensitive(sort_list):
    return sorted(sort_list, key=str.casefold)


def dir_check(path):
    if not exists(path):
        makedirs(path)


def run(home, providers, func, params, desc):
    for i in tqdm(range(0, len(providers)), desc=desc):
        provider = providers[i]
        func(home, provider, params)


def create_cols_in_df(df, columns_to_ensure):
    # Loop through the column names in columns_to_ensure
    for col_name in columns_to_ensure:
        if col_name not in df.columns:
            # Create an empty column with NaN values
            df[col_name] = ""  # You can also use df[col_name] = pd.Series(dtype='float64')
    return df

# https://ftp.ensembl.org/pub/current_json/homo_sapiens/homo_sapiens.json
# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/

def extract_ensembleid_from_dbxrefs(row):
    row = row.split("|")
    id = [x.split(":")[1] for x in row if "Ensembl" in x]
    if len(id)>0:
        return id[0]
    else:
        return ""

def extract_hgncid_from_dbxrefs(row):
    row = row.split("|")
    id = [x.split(":")[2] for x in row if "HGNC" in x]
    if len(id)>0:
        return "HGNC:"+str(id[0])
    else:
        return ""

def convert_cytoband2coord(row, cytobands):
    if not row.start > 0:
        sample = cytobands[cytobands.chromosome == "chr"+str(row.chromosome)]
        sample.cytoband = sample.chromosome.replace("chr", "", regex=True) + sample.cytoband
        if len(sample[sample.cytoband == row.map_location])==1:
                    row['start'] = sample[sample.cytoband == row.map_location]['start_pos'].reset_index(drop=True)[0]
                    row['end'] = sample[sample.cytoband == row.map_location]['end_pos'].reset_index(drop=True)[0]
                    row['strand'] = 1.0
    return row

def get_location_from_synonym(row):
    gene_map = {"BAT-25": "", "BAT-26": "", "CDKN2Ap14ARF": "", "CDKN2Ap16INK4A": "",
                "D17S250": "", "D2S123": "", "NR-21" : "", "NR-24": "",
                "LOC80054": ""}
    if row['symbol'] in gene_map.keys() and False:
        symbol = gene_map[row['symbol']]
    else:
        symbol = row['symbol']
    if isinstance(symbol, str):
        pattern = "(^|\|)"+symbol+"($|\|)"
        match = geneSymbol_location.loc[geneSymbol_location.Synonyms.str.contains(pattern)].reset_index(drop=True)
        if len(match) > 1:
            match = match.iloc[:1]
        if len(match) > 0:
            row['chromosome'], row['strand'], row['seq_start_position'], row['seq_end_position'], row['ncbi_gene_id'], row['ensembl_gene_id'] = match['chromosome'][0], match['strand'][0], match['start'][0], match['end'][0], match['GeneID'][0], match['ensembl_id'][0]
    return row

def get_geneSymbol_locations(paths):
    Reference = pd.read_csv(paths[0], sep='\t')
    #data = json.load( open(paths[0]))
    #data = data['genes']
    #Reference = pd.DataFrame(data)[['id', 'name', 'description', 'biotype', 'seq_region_name', 'strand', 'start', 'end', 'coord_system', 'HGNC', 'xrefs', 'RefSeq_mRNA']]
    #Reference = Reference[["id","name", "seq_region_name", "strand", "start", "end", "coord_system", "synonyms"]]
    #Reference = pd.concat([Reference, Reference.coord_system.apply(pd.Series)], axis=1).drop("coord_system", axis=1)
    #Reference["symbol"] = Reference.iloc[:,1]
    #Reference = Reference[Reference["symbol"].isna() == False]
    #Reference = Reference[Reference["version"] == "GRCh38"]

    NCBI_ref = pd.read_csv(paths[1], sep='\t')
    NCBI_ref['ensembl_id'] = NCBI_ref.dbXrefs.apply(lambda x: extract_ensembleid_from_dbxrefs(x))
    NCBI_ref['hgnc_id'] = NCBI_ref.dbXrefs.apply(lambda x: extract_hgncid_from_dbxrefs(x))
    NCBI_ref = NCBI_ref[["Symbol", "Synonyms", "chromosome", "map_location", "GeneID", "ensembl_id"]]

    cyto2coordinates = pd.read_csv(paths[2], sep='\t', names=["chromosome", "start_pos", "end_pos", "cytoband", "info"])
    cyto2coordinates = cyto2coordinates[cyto2coordinates.cytoband.isna() == False]

    GeneSymbol_Locations = NCBI_ref.merge(Reference, left_on="ensembl_id", right_on="id", how="left").apply(convert_cytoband2coord, cytobands=cyto2coordinates, axis=1)
    GeneSymbol_Locations = GeneSymbol_Locations[["Symbol", "Synonyms", "chromosome", "strand", "start", "end", "GeneID", "ensembl_id"]]
    GeneSymbol_Locations = GeneSymbol_Locations[GeneSymbol_Locations.start.isna() ==False]
    GeneSymbol_Locations = GeneSymbol_Locations.drop_duplicates(subset=['Symbol'])
    GeneSymbol_Locations.to_csv("../resources/gene_lcoation.tsv", sep='\t', index=False)

def merged_metadata(path):
    tsv_files = sorted([join(path, f) for f in get_files(path) if f.endswith('.tsv') and not f.__contains__('image') and not f.__contains__('molecular_metadata')])
    data = []
    for f in tsv_files:
        metadata = read_metadata_without_fields(f).fillna('Not Provided')
        metadata = metadata.applymap(lambda x: x.lower() if type(x) == str else str(int(x)))
        metadata = metadata.fillna('not provided').replace(
            ['unspecified', 'not reported', 'not specified', 'not collected', 'unknown'], 'not provided')
        if f.__contains__('patient.tsv'):
            metadata = metadata.replace('self reported', 'self-assessed')
            metadata['age_at_initial_diagnosis'] = [r if type(r) == str else str(int(r)) for r in metadata['age_at_initial_diagnosis']]
            data.append(metadata)
        elif f.__contains__('patient_sample.tsv'):
            data.append(metadata)
        elif f.__contains__('pdx_model.tsv'):
            data.append(metadata)
        elif f.__contains__('model_validation.tsv'):
            data.append(metadata)
        elif f.__contains__('cell_model.tsv'):
            data.append(metadata)
        elif f.__contains__('sharing.tsv'):
            data.append(metadata)
    return data


def generate_pdcm_sheet(df, path, provider, sheet_name):
    sheet_path = join(templates, 'metadata', f'metadata_template-{sheet_name}.tsv')
    sheet = read_metadata_with_fields(sheet_path)
    for col in sheet.columns:
        if col not in df.columns:
            df[col] = ''
    temp = pd.concat([sheet, df[sheet.columns]])
    temp.to_csv(join(path, f'{provider}_metadata-{sheet_name}.tsv'), sep='\t', index=False)

def input_to_pdcm(df, out_path, provider):
    df['Field'] = ''
    sheets = ['patient', 'patient_sample', 'model_validation', 'pdx_model', 'cell_model', 'sharing']
    for sheet in sheets:
        generate_pdcm_sheet(df, out_path, provider, sheet)