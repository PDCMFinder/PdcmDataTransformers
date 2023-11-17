from os import listdir, makedirs, getcwd
from os.path import isfile, join, isdir, exists
import pandas as pd
from tqdm import tqdm
from math import isnan

geneSymbol_location = pd.read_csv("/Users/tushar/CancerModels/utils/PdcmDataTransformers/resources/genes.tsv", sep='\t')

def get_dirs(path):
    return [f for f in listdir(path) if isdir(join(path, f))]


def get_files(path):
    return [f for f in listdir(path) if isfile(join(path, f))]


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
        else:
            print(symbol)
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
