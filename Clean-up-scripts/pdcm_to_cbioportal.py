# %%
from os import listdir, getcwd, rename, makedirs, remove
from os.path import isfile, join, isdir, exists
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
import subprocess as sp
from fuzzywuzzy import process, fuzz


# %% md
### Common functions
# %%
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
    hugo2ncbi = pd.read_csv("hugo_to_ncbi.txt", sep='\t')
    hugo2ncbi["NCBI Gene ID"] = hugo2ncbi["NCBI Gene ID"].fillna(0).astype(int)
    hugo2ncbi["NCBI Gene ID(supplied by NCBI)"] = hugo2ncbi["NCBI Gene ID(supplied by NCBI)"].fillna(0).astype(int)

    hugo2ncbi["NCBI Gene ID"] = np.where(hugo2ncbi["NCBI Gene ID"] != 0, hugo2ncbi["NCBI Gene ID"],
                                         hugo2ncbi["NCBI Gene ID(supplied by NCBI)"])
    hugo2ncbi['RefSeq IDs'] = np.where(hugo2ncbi["RefSeq IDs"].fillna('') != '', hugo2ncbi["RefSeq IDs"],
                                       hugo2ncbi["RefSeq(supplied by NCBI)"].fillna(''))

    hugo2ncbi = hugo2ncbi[hugo2ncbi["NCBI Gene ID"] != 0]
    return hugo2ncbi


# %% md
## Generate data for cBioPortal
# %% md
### Generate case lists
# %%
import os
import os.path
import sys

CASE_LIST_CONFIG_HEADER_COLUMNS = ["CASE_LIST_FILENAME", "STAGING_FILENAME", "META_STABLE_ID",
                                   "META_CASE_LIST_CATEGORY", "META_CANCER_STUDY_ID", "META_CASE_LIST_NAME",
                                   "META_CASE_LIST_DESCRIPTION"]
CASE_LIST_UNION_DELIMITER = "|"
CASE_LIST_INTERSECTION_DELIMITER = "&"
MUTATION_STAGING_GENERAL_PREFIX = "data_mutations"
SEQUENCED_SAMPLES_FILENAME = "sequenced_samples.txt"
MUTATION_CASE_LIST_META_HEADER = "sequenced_samples"
MUTATION_CASE_ID_COLUMN_HEADER = "Tumor_Sample_Barcode"
SAMPLE_ID_COLUMN_HEADER = "SAMPLE_ID"
NON_CASE_IDS = frozenset(
    ["MIRNA", "LOCUS", "ID", "GENE SYMBOL", "ENTREZ_GENE_ID", "HUGO_SYMBOL", "LOCUS ID", "CYTOBAND",
     "COMPOSITE.ELEMENT.REF", "HYBRIDIZATION REF"])
CANCER_STUDY_TAG = "<CANCER_STUDY>"
NUM_CASES_TAG = "<NUM_CASES>"


def generate_case_lists(case_list_config_filename, case_list_dir, study_dir, study_id, overwrite=False, verbose=False):
    header = []
    with open(case_list_config_filename, 'r') as case_list_config_file:
        # get header and validate
        header = case_list_config_file.readline().rstrip('\n').rstrip('\r').split('\t')
        # check full header matches what we expect
        for column in CASE_LIST_CONFIG_HEADER_COLUMNS:
            if column not in header:
                print >> sys.stderr, "ERROR: column '%s' is not in '%s'" % (column, case_list_config_filename)
                sys.exit(2)

        for line in case_list_config_file:
            line = line.rstrip('\n').rstrip('\r')
            config_fields = line.split('\t')
            case_list_filename = config_fields[header.index("CASE_LIST_FILENAME")]
            staging_filename_list = config_fields[header.index("STAGING_FILENAME")]
            case_list_file_full_path = os.path.join(case_list_dir, case_list_filename)
            if os.path.isfile(case_list_file_full_path) and not overwrite:
                if verbose:
                    print("LOG: generate_case_lists(), '%s' exists and overwrite is false, skipping caselist..." % (
                        case_list_filename))
                continue

            # might be single staging file
            staging_filenames = []
            # union (like all cases)
            union_case_list = CASE_LIST_UNION_DELIMITER in staging_filename_list
            # intersection (like complete or cna-seq)
            intersection_case_list = CASE_LIST_INTERSECTION_DELIMITER in staging_filename_list
            delimiter = CASE_LIST_UNION_DELIMITER if union_case_list else CASE_LIST_INTERSECTION_DELIMITER
            staging_filenames = staging_filename_list.split(delimiter)
            if verbose:
                print("LOG: generate_case_lists(), staging filenames: %s" % (",".join(staging_filenames)))

            # if this is intersection all staging files must exist
            if intersection_case_list and \
                    not all([os.path.isfile(os.path.join(study_dir, intersection_filename)) for intersection_filename in
                             staging_filenames]):
                continue

            # this is the set we will pass to write_case_list_file
            case_set = set([])
            # this indicates the number of staging files processed -
            # used to verify that an intersection should be written
            num_staging_files_processed = 0
            for staging_filename in staging_filenames:
                if verbose:
                    print("LOG: generate_case_lists(), processing staging file '%s'" % (staging_filename))
                # compute the case set
                case_list = []
                case_list = get_case_list_from_staging_file(study_dir, staging_filename, verbose)

                if len(case_list) == 0:
                    if verbose:
                        print("LOG: generate_case_lists(), no cases in '%s', skipping..." % (staging_filename))
                    continue

                if intersection_case_list:
                    if len(case_set) == 0:
                        # it is empty so initialize it
                        case_set = set(case_list)
                    else:
                        case_set = case_set.intersection(case_list)
                else:
                    # union of files or single file
                    case_set = case_set.union(case_list)

                num_staging_files_processed += 1

            # write case list file (don't make empty case lists)
            if len(case_set) > 0:
                if verbose:
                    print("LOG: generate_case_lists(), calling write_case_list_file()...")

                # do not write out complete cases file unless we've processed all the files required
                if intersection_case_list and num_staging_files_processed != len(staging_filenames):
                    if verbose:
                        print(
                            "LOG: generate_case_lists(), number of staging files processed (%d) != number of staging files required (%d) for '%s', skipping call to write_case_list_file()..." % (
                            num_staging_files_processed, len(staging_filenames), case_list_filename))
                else:
                    write_case_list_file(header, config_fields, study_id, case_list_file_full_path, case_set, verbose)
            elif verbose:
                print("LOG: generate_case_lists(), case_set.size() == 0, skipping call to write_case_list_file()...")


def get_case_list_from_staging_file(study_dir, staging_filename, verbose):
    if verbose:
        print("LOG: get_case_list_from_staging_file(), '%s'" % (staging_filename))

    case_set = set([])

    # if we are processing mutations data and a SEQUENCED_SAMPLES_FILENAME exists, use it
    if MUTATION_STAGING_GENERAL_PREFIX in staging_filename:
        sequenced_samples_full_path = os.path.join(study_dir, SEQUENCED_SAMPLES_FILENAME)
        if os.path.isfile(sequenced_samples_full_path):
            if verbose:
                print(
                    "LOG: get_case_list_from_staging_file(), '%s' exists, calling get_case_list_from_sequenced_samples_file()" % (
                        SEQUENCED_SAMPLES_FILENAME))
            return get_case_list_from_sequenced_samples_file(sequenced_samples_full_path, verbose)

    staging_file_full_path = os.path.join(study_dir, staging_filename)
    if not os.path.isfile(staging_file_full_path):
        return []

    # staging file
    with open(staging_file_full_path, 'r') as staging_file:
        id_column_index = 0
        process_header = True
        for line in staging_file:
            line = line.rstrip('\n')
            if line.startswith('#'):
                if line.startswith('#' + MUTATION_CASE_LIST_META_HEADER + ':'):
                    # split will split on any whitespace, tabs or any number of consecutive spaces
                    return line[len(MUTATION_CASE_LIST_META_HEADER) + 2:].strip().split()
                continue  # this is a comment line, skip it
            values = line.split('\t')

            # is this the header line?
            if process_header:
                # look for MAF file case id column header
                # if this is not a MAF file and header contains the case ids, return here
                # we are assuming the header contains the case ids because SAMPLE_ID_COLUMN_HEADER is missing
                if MUTATION_CASE_ID_COLUMN_HEADER not in values and SAMPLE_ID_COLUMN_HEADER not in [x.upper() for x in
                                                                                                    values]:
                    if verbose:
                        print(
                            "LOG: get_case_list_from_staging_file(), this is not a MAF header but has no '%s' column, we assume it contains sample ids..." % (
                                SAMPLE_ID_COLUMN_HEADER))
                    for potential_case_id in values:
                        # check to filter out column headers other than sample ids
                        if potential_case_id.upper() in NON_CASE_IDS:
                            continue
                        case_set.add(potential_case_id)
                    break  # got case ids from header, don't read the rest of the file
                else:
                    # we know at this point one of these columns exists, so no fear of ValueError from index method
                    id_column_index = values.index(
                        MUTATION_CASE_ID_COLUMN_HEADER) if MUTATION_CASE_ID_COLUMN_HEADER in values else [x.upper() for
                                                                                                          x in
                                                                                                          values].index(
                        SAMPLE_ID_COLUMN_HEADER)
                    if verbose:
                        print(
                            "LOG: get_case_list_from_staging_file(), this is a MAF or clinical file, samples ids in column with index: %d" % (
                                id_column_index))
                process_header = False
                continue  # done with header, move on to next line
            case_set.add(values[id_column_index])

    return list(case_set)


def get_case_list_from_sequenced_samples_file(sequenced_samples_full_path, verbose):
    if verbose:
        print("LOG: get_case_list_from_sequenced_samples_file, '%s'", sequenced_samples_full_path)

    case_set = set([])
    with open(sequenced_samples_full_path, 'r') as sequenced_samples_file:
        for line in sequenced_samples_file:
            case_set.add(line.rstrip('\n'))

    if verbose:
        print("LOG: get_case_list_from_sequenced_samples_file, case set size: %d" % (len(case_set)))

    return list(case_set)


def write_case_list_file(case_list_config_header, case_list_config_fields, study_id, case_list_full_path, case_set,
                         verbose):
    if verbose:
        print("LOG: write_case_list_file(), '%s'" % (case_list_full_path))
    with open(case_list_full_path, 'w') as case_list_file:
        case_list_file.write("cancer_study_identifier: " + study_id + "\n")
        stable_id = case_list_config_fields[case_list_config_header.index("META_STABLE_ID")].replace(CANCER_STUDY_TAG,
                                                                                                     study_id)
        case_list_file.write("stable_id: " + stable_id + "\n")
        case_list_file.write(
            "case_list_name: " + case_list_config_fields[case_list_config_header.index("META_CASE_LIST_NAME")] + "\n")
        case_list_description = case_list_config_fields[
            case_list_config_header.index("META_CASE_LIST_DESCRIPTION")].replace(NUM_CASES_TAG, str(len(case_set)))
        case_list_file.write("case_list_description: " + case_list_description + "\n")
        case_list_file.write("case_list_category: " + case_list_config_fields[
            case_list_config_header.index("META_CASE_LIST_CATEGORY")] + "\n")
        case_list_file.write("case_list_ids: " + '\t'.join(case_set) + "\n")


def main(case_list_config_filename, case_list_dir, study_dir, study_id, overwrite, verbose):
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--case-list-config-file', action = 'store', dest = 'case_list_config_file', required = True, help = 'Path to the case list configuration file.  An example can be found in "test/resources/generate_case_lists/case_list_config.tsv"')
    parser.add_argument('-d', '--case-list-dir', action = 'store', dest = 'case_list_dir', required = True, help = 'Path to the directory in which the case list files should be written')
    parser.add_argument('-s', '--study-dir', action = 'store', dest = 'study_dir', required = True, help = 'The directory that contains the cancer study genomic files')
    parser.add_argument('-i', '--study-id', action = 'store', dest = 'study_id', required = True, help = 'The cancer study stable id')
    parser.add_argument('-o', '--overwrite', action = 'store_true', dest = 'overwrite', required = False, help = 'When given, overwrite the case list files')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', required = False, help = 'When given, be verbose')
    args = parser.parse_args()

    #case_list_config_filename = args.case_list_config_file
    #case_list_dir = args.case_list_dir
    #study_dir = args.study_dir
    #study_id = args.study_id
    overwrite = args.overwrite
    verbose = args.verbose
    '''
    if verbose:
        print("LOG: case_list_config_file='%s'" % (case_list_config_filename))
        print("LOG: case_list_dir='%s'" % (case_list_dir))
        print("LOG: study_dir='%s'" % (study_dir))
        print("LOG: study_id='%s'" % (study_id))
        print("LOG: overwrite='%s'" % (overwrite))
        print("LOG: verbose='%s'" % (verbose))

    if not os.path.isfile(case_list_config_filename):
        print("ERROR: case list configuration file '%s' does not exist or is not a file" % (case_list_config_filename),
              file=sys.stderr)
        sys.exit(2)

    if not os.path.isdir(case_list_dir):
        print("ERROR: case list file directory '%s' does not exist or is not a directory" % (case_list_dir),
              file=sys.stderr)
        sys.exit(2)

    if not os.path.isdir(study_dir):
        print("ERROR: study directory '%s' does not exist or is not a directory" % (study_dir), file=sys.stderr)
        sys.exit(2)

    generate_case_lists(case_list_config_filename, case_list_dir, study_dir, study_id, overwrite, verbose)


# %% md
### Clinical files
# %%
def add_headers_clinical_patient():
    # Column headers - The attribute Display Names: The display name for each clinical attribute
    column_headers = ["#Patient Identifier", "Sex", "Diagnosis Age", "Overall Survival (Months)",
                      "Overall Survival Status"]
    # Column descriptions - The attribute Descriptions: Long(er) description of each clinical attribute
    column_description = ["#Identifier to uniquely specify a patient.", "Sex",
                          "Age at which a condition or disease was first diagnosed.",
                          "Overall survival in months since initial diagonosis.", "Overall patient survival status."]
    # Column data type - The attribute Datatype: The datatype of each clinical attribute (must be one of: STRING, NUMBER, BOOLEAN)
    column_data_type = ["#STRING", "STRING", "NUMBER", "NUMBER", "STRING"]
    # Column priority - A number which indicates the importance of each attribute. A higher number indicates a higher priority.
    column_priority = ["#1", "1", "1", "1", "1"]
    # Column headers for validation
    column_header_validation = ["PATIENT_ID", "SEX", "AGE", "OS_MONTHS", "OS_STATUS"]
    return [column_headers, column_description, column_data_type, column_priority, column_header_validation]


def add_headers_clinical_sample():
    # Column headers - The attribute Display Names: The display name for each clinical attribute
    column_headers = ["#Patient Identifier", "Sample Identifier", "Tumor Type", "Cancer Type", "Cancer Type Detailed",
                      "Primary site", "Tumor Grade", "Model type", "Model ID"]
    # Column descriptions - The attribute Descriptions: Long(er) description of each clinical attribute
    column_description = ["#Identifier to uniquely specify a patient.", "A unique sample identifier.",
                          "The type of tumour sample (i.e., normal, primary, met, recurrence).", "Cancer Type",
                          "Cancer Type Detailed",
                          "Site of the primary tumor where primary cancer is originating from (may not correspond to the site of the current tissue sample). No abbreviations.",
                          "The implanded tumour grade value.", "Type of patient derived cancer model",
                          "Unique identifier for the PDCMs"]
    # Column data type - The attribute Datatype: The datatype of each clinical attribute (must be one of: STRING, NUMBER, BOOLEAN)
    column_data_type = ["#STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING"]
    # Column priority - A number which indicates the importance of each attribute. A higher number indicates a higher priority.
    column_priority = ["#1", "1", "1", "1", "1", "1", "1", "1", "1"]
    # Column headers for validation
    column_header_validation = ["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE", "CANCER_TYPE", "CANCER_TYPE_DETAILED",
                                "PRIMARY_SITE", "TUMOR_GRADE", "MODEL_TYPE", "MODEL_ID"]
    return [column_headers, column_description, column_data_type, column_priority, column_header_validation]


def generate_clinical_df(headers):
    # Generate the df
    c_bio_data_clinical = pd.DataFrame(columns=headers[0])
    c_bio_data_clinical.loc[0, :] = headers[1]
    c_bio_data_clinical.loc[1, :] = headers[2]
    c_bio_data_clinical.loc[2, :] = headers[3]
    c_bio_data_clinical.loc[3, :] = headers[4]
    return c_bio_data_clinical


def write_clinical_file(out_df, out_path, file_type):
    out_file = join(out_path, "data_clinical_" + file_type + ".txt")
    out_df.to_csv(out_file, sep="\t", index=False)


def clinical_patient_data_transformation(data_path, provider, out_path):
    out_df = generate_clinical_df(add_headers_clinical_patient())
    data_path = join(data_path, provider + "_metadata-patient.tsv")

    mapper = {"patient_id": "#Patient Identifier", "sex": "Sex", "age_at_initial_diagnosis": "Diagnosis Age"}
    temp = read_metadata_without_fields(data_path).rename(columns=mapper)

    temp["Overall Survival (Months)"] = ""
    temp["Overall Survival Status"] = ""
    cs = pd.read_csv(join(out_path, "data_clinical_sample.txt"), sep="\t")
    patient_ids = list(cs["#Patient Identifier"])
    temp = temp[temp["#Patient Identifier"].isin(patient_ids)]
    # temp = temp.merge(cs[["#Patient Identifier", "Cancer Type"]], on="#Patient Identifier", how="inner")
    # temp["Tumor Type"] = temp["Cancer Type"]
    out_df = pd.concat([out_df, temp[out_df.columns]], ignore_index=True)

    return out_df.replace("(?i)not provided", "", regex=True).replace("(?i)not collected", "",
                                                                      regex=True).drop_duplicates(
        subset=["#Patient Identifier"])


def clinical_patient_sample_data_transformation(data_path, provider):
    out_df = generate_clinical_df(add_headers_clinical_sample())
    data_path = join(data_path, provider + "_metadata-patient_sample.tsv")

    mapper = {"patient_id": "#Patient Identifier", "sample_id": "Sample Identifier", "tumour_type": "Tumor Type OG",
              "primary_site": "Primary site OG", "grade": "Tumor Grade"}
    temp = read_metadata_without_fields(data_path).rename(columns=mapper)
    temp = temp.merge(get_meta_from_api(provider), on="model_id", how="inner").reset_index(drop=True)
    temp["Cancer Type"] = temp["cancer_system"]
    temp["Cancer Type Detailed"] = temp["histology"]
    temp["Model type"] = temp["type"]
    temp["Tumor Type"] = temp["tumor_type"]
    temp["Primary site"] = temp["primary_site"]
    temp["Model ID"] = temp["model_id"]
    temp["Tumor Grade"] = temp["Tumor Grade"].str.replace('Not provided', '', flags=re.IGNORECASE).str.replace(
        'Not collected', '', flags=re.IGNORECASE).replace("p.*", "", regex=True).replace("T.*", "", regex=True).replace(
        " \(.*", "", regex=True).replace(";.*", "", regex=True).str.replace('Not collected', '',
                                                                            flags=re.IGNORECASE).replace("Moderate",
                                                                                                         "moderately differentiated",
                                                                                                         regex=True).replace(
        "Poorly", "poorly differentiated", regex=True).replace("Well", "well differentiated", regex=True).replace(
        "High-Grade", "high", regex=True).replace("High grade", "high", regex=True).replace("Grade ", "G", regex=True)
    out_df = pd.concat([out_df, temp[out_df.columns]], ignore_index=True)
    return out_df.replace("(?i)not provided", "", regex=True).replace("(?i)not collected", "",
                                                                      regex=True).drop_duplicates(
        subset=["#Patient Identifier", "Sample Identifier"])


def get_meta_from_api(provider_name):
    df = pd.read_json("https://www.cancermodels.org/api/model_metadata?data_source=eq." + provider_name)
    return df[["model_id", "cancer_system", "histology", "type", "tumor_type", "primary_site"]]


def generate_clinical_patient(in_path, out_path, provider):
    write_clinical_file(clinical_patient_data_transformation(in_path, provider, out_path), out_path, "patient")


def generate_clinical_sample(in_path, out_path, provider):
    write_clinical_file(clinical_patient_sample_data_transformation(in_path, provider), out_path, "sample")


# %% md
### Molecular data
# %%
def get_platform(path, type):
    df = read_metadata_without_fields(path).fillna("")
    df = df[df["molecular_characterisation_type"] == type]
    df["platform_name"] = df["library_strategy"] + " -  " + df["instrument_model"]
    if path.__contains__("JAX"):
        if type == "mutation":
            return "WES"
        if type == "expression":
            return "RNA-Seq"
        if type == "copy number alteration":
            return "SNP"
    return ", ".join(df["platform_name"])


def filter_sample_ids(x, sample_ids):
    out = [val for val in sample_ids if val in x]
    if len(out) > 0 and out[0] != '':
        return out[0]
    else:
        "No match"


def filter_samples(df, col, out_path):
    sample_ids = list(pd.read_csv(join(out_path, "data_clinical_sample.txt"), sep="\t")["Sample Identifier"])[4:]
    df = df[df[col].apply(lambda x: any(val in x for val in sample_ids))]
    df[col] = df[col].apply(lambda x: filter_sample_ids(x, sample_ids))
    df = df[df[col] != "No match"]
    df = df[df[col].fillna('') != ""]
    return df


def read_mol_data(path):
    files = get_files(path)
    if len(files) == 0:
        dirs = get_dirs(path)
        for dir in dirs:
            files.append(get_files(join(path, dir)))
        files = [item for sublist in files for item in sublist]
    df = pd.DataFrame()
    for file in files:
        temp = pd.read_csv(file, sep="\t")
        df = pd.concat([df, temp], ignore_index=True)
    return df


def compute_end_pos(df):
    df["End_Position"] = pd.to_numeric(df["Start_Position"]) + np.where(
        df.Reference_Allele.str.len() >= df.Tumor_Seq_Allele2.str.len(), df.Reference_Allele.str.len(),
        df.Tumor_Seq_Allele2.str.len()) - 1
    return df


def map_gene_symbol_to_id(symbol):
    # Check if the symbol is in the reference DataFrame
    if symbol in reference_df['Approved symbol'].tolist():
        return reference_df.loc[reference_df['Approved symbol'] == symbol, 'NCBI Gene ID'].values[0]
    # Check if the symbol is in the 'previous' column
    elif reference_df['Previous symbols'].str.contains(symbol, na=False).any():
        return reference_df.loc[reference_df['Previous symbols'].str.contains(symbol, na=False), 'NCBI Gene ID'].values[
            0]
    # Check if the symbol is in the 'alias' column
    elif reference_df['Alias symbols'].str.contains(symbol, na=False).any():
        return reference_df.loc[reference_df['Alias symbols'].str.contains(symbol, na=False), 'NCBI Gene ID'].values[0]
    else:
        return ""


def map_id_to_hugo(id, symbol):
    if id == "":
        return symbol
    return reference_df.loc[reference_df['NCBI Gene ID'] == id, 'Approved symbol'].values[0]


def map_id_to_refseq(id):
    if id == "":
        return ""
    else:
        return reference_df.loc[reference_df['NCBI Gene ID'] == id, 'RefSeq IDs'].values[0]


def handle_frameshift(row):
    if not row["Consequence"].__contains__("frameshift"):
        return row
    variant = row["variant_class"]
    if variant.__contains__("frameshift"):
        if row["coding_sequence_change"].__contains__("ins"):
            variant = "insertion"
        elif row["coding_sequence_change"].__contains__("del"):
            variant = "deletion"
    row["Consequence"] = "frameshift_variation_" + variant
    return row


def generate_variation_classification(df):
    mapper = {"frameshift_variant_deletion": "Frame_Shift_Del",
              "frameshift_variant_insertion": "Frame_Shift_Ins",
              "inframe_deletion": "In_Frame_Del",
              "inframe_insertion": "In_Frame_Ins",
              "missense_variant": "Missense_Mutation",
              "stop_gained": "Nonsense_Mutation",
              "5_prime_UTR_variant": "5'UTR",
              "upstream_gene_variant": "5'Flank",
              "downstream_gene_variant": "3'Flank",
              "3_prime_UTR_variant": "3'UTR",
              "non_coding_transcript_variant": "RNA",
              "intron_variant": "Intron",
              "splice_region_variant": "Splice_Region",
              "synonymous_variant": "Silent",
              "stop_lost": "Nonstop_Mutation",
              "start_retained_variant": "Translation_Start_Site",
              "intergenic_variant": "IGR"}
    # Targeted_Region, De_novo_Start_InFrame, De_novo_Start_OutOfFrame, Splice_Site (dup) }
    # df["temp_cons"] = np.where(df["Consequence"].str.contains("frameshift"), df["Consequence"], df["Consequence"])
    df = df.apply(lambda x: handle_frameshift(x), axis=1)
    df["Variant_Classification"] = df["Consequence"].apply(lambda x: mapper.get(x, "Unknown"))
    # np.where(df["Consequence"].to_list() in mapper.keys(), df["Consequence"], "Unknown"))
    df["Variant_Classification"] = df["Variant_Classification"].replace(mapper)
    return df


def generate_meta_mutation(study, out_path, platform):
    meta_df = pd.DataFrame(columns=[0, 1])
    meta_df.loc[0, :] = ["cancer_study_identifier:", study]
    meta_df.loc[1, :] = ["genetic_alteration_type:", "MUTATION_EXTENDED"]
    meta_df.loc[2, :] = ["datatype:", "MAF"]
    meta_df.loc[3, :] = ["stable_id:", "mutations"]
    meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "true"]
    meta_df.loc[5, :] = ["profile_description:", "Mutation data from " + platform]
    meta_df.loc[6, :] = ["profile_name:", "Mutations"]
    meta_df.loc[7, :] = ["variant_classification_filter:", "De_novo_Start_InFrame, De_novo_Start_OutOfFrame, Unknown"]
    meta_df.loc[8, :] = ["data_filename:", "data_mutations_extended.txt"]
    meta_df.to_csv(join(out_path, "meta_mutations_extended.txt"), sep="\t", index=False, header=False)


def generate_mutation_data(path, out_path):
    mapper = {"sample_id": "Tumor_Sample_Barcode", "seq_start_position": "Start_Position", "chromosome": "Chromosome",
              "ref_allele": "Reference_Allele", "alt_allele": "Tumor_Seq_Allele2", "symbol": "Hugo_Symbol",
              "consequence": "Consequence", "amino_acid_change": "HGVSp_Short",
              "ensembl_gene_id": "Gene", "ensembl_transcript_id": "Transcript_ID", "codon_change": "Codons",
              "ncbi_transcript_id": "RefSeq"}
    mut_df = read_mol_data(path).rename(columns=mapper).fillna("")
    mut_df = filter_samples(mut_df, "Tumor_Sample_Barcode", out_path)
    out_column = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position",
                  "End_Position", "Strand", "Consequence", "Variant_Classification", "Variant_Type", "Reference_Allele",
                  "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode",
                  "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
                  "Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
                  "Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status",
                  "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", "BAM_File", "Sequencer",
                  "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count", "HGVSc", "HGVSp", "HGVSp_Short",
                  "Transcript_ID", "RefSeq", "Protein_position", "Codons", "Hotspot"]
    mut_df = compute_end_pos(mut_df)
    mut_df["NCBI_Build"] = "GRCh38"
    mut_df["Strand"] = "+"
    mut_df["HGVSc"] = mut_df["Transcript_ID"] + "c." + mut_df["coding_sequence_change"]
    mut_df['Protein_position'] = np.where(mut_df['HGVSp_Short'].fillna('') != '', mut_df['HGVSp_Short'].str[1:-1], '')
    mut_df["HGVSp_Short"] = np.where(mut_df["HGVSp_Short"].fillna('') != "",
                                     "p." + mut_df["HGVSp_Short"].str.replace('.0', ''), "")
    mut_df["Entrez_Gene_Id"] = mut_df['Hugo_Symbol'].apply(map_gene_symbol_to_id)
    mut_df = mut_df[mut_df["Entrez_Gene_Id"] != ""]
    mut_df['Hugo_Symbol'] = mut_df[["Entrez_Gene_Id", 'Hugo_Symbol']].apply(lambda x: map_id_to_hugo(x[0], x[1]),
                                                                            axis=1)
    mut_df['RefSeq'] = mut_df['Entrez_Gene_Id'].apply(map_id_to_refseq)
    # mut_df["User_Amino_Acid_Change"] = mut_df["HGVSp_Short"]
    mut_df = generate_variation_classification(mut_df)

    mut_df["Variant_Type"] = mut_df["variant_class"].replace("SNV", "SNP").replace("insertion", "INS").replace(
        "deletion", "DEL")
    for column in out_column:
        if column not in mut_df.columns:
            mut_df[column] = ""
    if mut_df.shape[0] > 0:
        mut_df[out_column].to_csv(join(out_path, "data_mutations_extended.txt"), index=False, sep="\t")
    else:
        remove(join(out_path, "meta_mutations_extended.txt"))


def generate_meta_cna(study, out_path, platform, datatype):
    meta_df = pd.DataFrame(columns=[0, 1])
    meta_df.loc[0, :] = ["cancer_study_identifier:", study]
    meta_df.loc[1, :] = ["genetic_alteration_type:", "COPY_NUMBER_ALTERATION"]
    if datatype == "log":
        meta_df.loc[2, :] = ["datatype:", "LOG2-VALUE"]
        meta_df.loc[3, :] = ["stable_id:", "log2CNA"]
        out_name = "log2_cna"
        meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "false"]
        meta_df.loc[5, :] = ["profile_description:",
                             "Putative copy-number from GISTIC 2.0. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification. " + platform]
        meta_df.loc[6, :] = ["profile_name:", "Putative copy-number alterations from GISTIC"]
        meta_df.loc[7, :] = ["data_filename:", "data_log2_cna.txt"]

    elif datatype == "gistic":
        meta_df.loc[2, :] = ["datatype:", "DISCRETE"]
        meta_df.loc[3, :] = ["stable_id:", "gistic"]
        out_name = "cna"
        meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "false"]
        meta_df.loc[5, :] = ["profile_description:", "Log2 copy-number data from " + platform]
        meta_df.loc[6, :] = ["profile_name:", "Log2 copy-number values"]
        meta_df.loc[7, :] = ["data_filename:", "data_cna.txt"]
    meta_df.to_csv(join(out_path, "meta_" + out_name + ".txt"), sep="\t", index=False, header=False)


def generate_cna_data(path, out_path, study, platform):
    df = read_mol_data(path)[["sample_id", "symbol", "log2r_cna", "gistic_value"]]
    col_value = "log2r_cna"
    datatype = "log"
    out_name = "log2_cna"
    df = filter_samples(df, "sample_id", out_path)
    df['Hugo_Symbol'] = df['symbol']
    if df[col_value].isna().all():
        col_value = "gistic_value"
        datatype = "gistic"
        out_name = "cna"
        df['gistic_value'] = df['gistic_value'].astype(int)
    elif df[col_value].isna().all():
        return None
    if datatype == "gistic":
        df = df.pivot_table(index='Hugo_Symbol', columns='sample_id', values=col_value, aggfunc='first').fillna(
            100).astype(int).replace(100, 'NA').reset_index()
    else:
        df = df.pivot_table(index='Hugo_Symbol', columns='sample_id', values=col_value, aggfunc='first').fillna(
            "NA").reset_index()
    out_cols = list(df.columns)
    out_cols.insert(1, 'Entrez_Gene_Id')
    df['Entrez_Gene_Id'] = df['Hugo_Symbol'].apply(map_gene_symbol_to_id)
    df = df[df["Entrez_Gene_Id"] != ""]
    df['Hugo_Symbol'] = df[["Entrez_Gene_Id", 'Hugo_Symbol']].apply(lambda x: map_id_to_hugo(x[0], x[1]), axis=1)
    if df.shape[0] > 0:
        generate_meta_cna(study, out_path, platform, datatype)
        df[out_cols].to_csv(join(out_path, "data_" + out_name + ".txt"), index=False, sep="\t")


def generate_meta_expression(study, out_path, platform, datatype):
    meta_df = pd.DataFrame(columns=[0, 1])
    meta_df.loc[0, :] = ["cancer_study_identifier:", study]
    meta_df.loc[1, :] = ["genetic_alteration_type:", "MRNA_EXPRESSION"]
    if datatype == "mrna":
        meta_df.loc[2, :] = ["datatype:", "CONTINUOUS"]
        meta_df.loc[3, :] = ["stable_id:", "rna_seq_mrna"]
    elif datatype == "Zscore":
        meta_df.loc[2, :] = ["datatype:", "Z-SCORE"]
        meta_df.loc[3, :] = ["stable_id:", "mrna_median_Zscores"]
    meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "true"]
    meta_df.loc[5, :] = ["profile_description:", "Expression data from " + platform]
    meta_df.loc[6, :] = ["profile_name:", "Expression"]
    meta_df.loc[7, :] = ["data_filename:", "data_mrna_seq_rpkm.txt"]
    meta_df.to_csv(join(out_path, "meta_mrna_seq_rpkm.txt"), sep="\t", index=False, header=False)


def generate_expression_data(path, out_path, platform, study):
    df = read_mol_data(path)[["sample_id", "symbol", "rnaseq_fpkm", "z_score"]]
    df = filter_samples(df, "sample_id", out_path)

    df["Hugo_Symbol"] = df["symbol"]
    no_fpkm = df['rnaseq_fpkm'].isna().all()
    datatype = "mrna"
    value_column = "rnaseq_fpkm"
    if no_fpkm:
        print("Using Z score: " + path)
        datatype = "Zscore"
        value_column = "z_score"
    df = df[["sample_id", "Hugo_Symbol", value_column]]
    df = df.pivot_table(index='Hugo_Symbol', columns='sample_id', values=value_column, aggfunc='first').fillna(
        'NA').reset_index()
    out_cols = list(df.columns)
    out_cols.insert(1, 'Entrez_Gene_Id')
    df['Entrez_Gene_Id'] = df['Hugo_Symbol'].apply(map_gene_symbol_to_id)
    df = df[df["Entrez_Gene_Id"] != ""]
    df['Hugo_Symbol'] = df[["Entrez_Gene_Id", 'Hugo_Symbol']].apply(lambda x: map_id_to_hugo(x[0], x[1]), axis=1)

    if df.shape[0] > 0:
        generate_meta_expression(study, out_path, platform, datatype)
        df[out_cols].to_csv(join(out_path, "data_mrna_seq_rpkm.txt"), sep="\t", index=False)


def generate_cna_files(in_path, study, provider, out_path):
    platform = get_platform(join(in_path, provider + "_molecular_metadata-platform.tsv"), "mutation")
    # generate_meta_cna(study, out_path, platform)
    generate_cna_data(join(in_path, "cna"), out_path, study, platform)


def generate_mutation_files(in_path, study, provider, out_path):
    platform = get_platform(join(in_path, provider + "_molecular_metadata-platform.tsv"), "mutation")
    generate_meta_mutation(study, out_path, platform)
    generate_mutation_data(join(in_path, "mut"), out_path)


def generate_expression_files(in_path, study, provider, out_path):
    platform = get_platform(join(in_path, provider + "_molecular_metadata-platform.tsv"), "expression")
    generate_expression_data(join(in_path, "expression"), out_path, platform, study)


# %% md
### Timeline data
# %%
def generate_meta_timeline(study, type, out_path):
    meta_df = pd.DataFrame(columns=[0, 1])
    meta_df.loc[0, :] = ["cancer_study_identifier:", study]
    meta_df.loc[1, :] = ["genetic_alteration_type:", "CLINICAL"]
    meta_df.loc[2, :] = ["datatype:", "TIMELINE"]
    meta_df.loc[7, :] = ["data_filename:", "data_timeline_" + type + ".txt"]
    meta_df.to_csv(join(out_path, "meta_timeline_" + type + ".txt"), sep="\t", index=False, header=False)


def generate_timeline_data_file(in_path, provider, study, out_path, type):
    ps = read_metadata_without_fields(join(in_path, provider + "_metadata-patient_sample.tsv"))
    if type == 'specimen':
        out_cols = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SAMPLE_ID", "SPECIMEN_SITE", "SOURCE"]
        ps["PATIENT_ID"], ps["SAMPLE_ID"], ps["STOP_DATE"] = ps["patient_id"], ps["sample_id"], ""
        ps['START_DATE'], ps['SOURCE'], ps['EVENT_TYPE'] = 0, provider, type.upper()
        ps["SPECIMEN_SITE"] = ps["collection_site"]
        ids = list(pd.read_csv(join(out_path, "data_clinical_patient.txt"), sep="\t")["#Patient Identifier"])
        ps = ps[ps["patient_id"].isin(ids)]
        ps[out_cols].to_csv(join(out_path, "data_timeline_" + type + ".txt"), sep="\t", index=False)
    if type == "treatment":
        treatment = read_metadata_with_fields(join(in_path, "treatment", provider + "_patienttreatment-Sheet1.tsv"))
        dates = ps[["patient_id", "collection_date"]]
        try:
            dates["date_time"] = pd.to_datetime(
                dates["collection_date"].replace("(?i)not provided", "", regex=True).replace("(?i)not collected", "",
                                                                                             regex=True).str.replace(
                    " ", "-"), format='%b-%y')
        except:
            try:
                dates["date_time"] = pd.to_datetime(
                    dates["collection_date"].replace("(?i)not provided", "", regex=True).replace("(?i)not collected",
                                                                                                 "",
                                                                                                 regex=True).str.replace(
                        " ", "-"), format='%b-%Y')
            except:
                try:
                    dates["date_time"] = pd.to_datetime(
                        dates["collection_date"].replace("(?i)not provided", "", regex=True).replace(
                            "(?i)not collected", "", regex=True).str.replace(" ", "-"), format='%d/%b/%y')
                except:
                    print("Failed to handle the dates")
                    dates["date_time"] = 0

        treatment = treatment.merge(dates, on="patient_id", how="inner").reset_index(drop=True)
        try:
            treatment["start_date_time"] = pd.to_datetime(
                treatment["treatment_starting_date"].str.replace('Not provided', '', flags=re.IGNORECASE).str.replace(
                    'Not collected', '', flags=re.IGNORECASE).str.replace(" ", "-"), format='%b-%y')
        except:
            try:
                treatment["start_date_time"] = pd.to_datetime(
                    treatment["treatment_starting_date"].str.replace('Not provided', '',
                                                                     flags=re.IGNORECASE).str.replace('Not collected',
                                                                                                      '',
                                                                                                      flags=re.IGNORECASE).str.replace(
                        " ", "-"), format='%b-%Y')
            except:
                try:
                    treatment["start_date_time"] = pd.to_datetime(
                        treatment["treatment_starting_date"].str.replace('Not provided', '',
                                                                         flags=re.IGNORECASE).str.replace(
                            'Not collected', '', flags=re.IGNORECASE).str.replace(" ", "-"), format='%d/%b/%y')
                except:
                    treatment["start_date_time"] = 0
        out_cols = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TREATMENT_TYPE", "AGENT"]
        treatment["PATIENT_ID"] = treatment["patient_id"]
        treatment["START_DATE"] = (treatment['start_date_time'] - treatment['date_time']).dt.days
        treatment["STOP_DATE"] = pd.to_numeric(treatment["START_DATE"], errors='coerce') + (
                    pd.to_numeric(treatment["treatment_duration"].replace("(?i)not provided", "0", regex=True),
                                  errors='coerce') * 30)
        treatment["START_DATE"].fillna(0, inplace=True)
        treatment["STOP_DATE"].fillna(0, inplace=True)
        treatment["EVENT_TYPE"] = type.upper()

        treatment["START_DATE"] = treatment["START_DATE"].astype(int)
        treatment["STOP_DATE"] = treatment["STOP_DATE"].astype(int)
        treatment["TREATMENT_TYPE"] = "Medical Therapy"
        treatment["AGENT"] = treatment["treatment_name"]

        treatment[out_cols].to_csv(join(out_path, "data_timeline_" + type + ".txt"), sep="\t", index=False)
    if type == "lab_test":
        cyto = read_metadata_with_fields(join(in_path, "cyto", provider + "_cytogenetics-Sheet1.tsv"))
        cyto = cyto.merge(ps[["patient_id", "sample_id"]], on="sample_id", how="inner")
        out_cols = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT"]
        cyto["PATIENT_ID"], cyto["START_DATE"], cyto["STOP_DATE"] = cyto["patient_id"], 0, ""
        cyto["EVENT_TYPE"] = type.upper()
        cyto["TEST"], cyto["RESULT"] = cyto['symbol'], cyto['marker_status']
        cyto[out_cols].to_csv(join(out_path, "data_timeline_" + type + ".txt"), sep="\t", index=False)


def generate_timeline_data(in_path, out_path, study, provider):
    generate_meta_timeline(study, "specimen", out_path)
    generate_timeline_data_file(in_path, provider, study, out_path, 'specimen')
    if exists(join(in_path, "treatment")):
        generate_meta_timeline(study, "treatment", out_path)
        generate_timeline_data_file(in_path, provider, study, out_path, 'treatment')
    if exists(join(in_path, "cyto")):
        generate_meta_timeline(study, "lab_test", out_path)
        generate_timeline_data_file(in_path, provider, study, out_path, 'lab_test')


# %% md
### Meta files
# %%
def get_provider_description(provider_name):
    if provider_name == "DFCI-CPDM":
        return "The Center for Patient Derived Models (CPDM) at Dana-Farber  Cancer Institute (DFCI) is a strategic collaborative research center with the expertise  to generate and characterize patient derived xenografts (PDX), patient derived cell  lines (PDCL - 3D organoid/spheroid and 2D adherent cultures), and acute cell models  drug testing. Through collaboration with major disease groups in the Dana-Farber Cancer Institute, Brigham and Women's Hospital, and Boston Children Hospital Cancer Centers, we have made a large collection of patient derived models of brain tumors, hematologic tumors, and many other solid tumors available to academic and industrial researchers worldwide."
    df = pd.read_json("https://www.cancermodels.org/api/provider_group?abbreviation=eq." + provider_name)
    description = str(df["description"][0]).replace("\n\n", " ").replace("\n", " ")
    if len(description) >= 1024:
        description = description[0:1020]
    return description


def get_study_cancer_type(provider):
    cancer_types = {"CCIA": "lymph", "CHOP": "nbl", "Curie-BC": "brca", "Curie-LC": "lung", "Curie-OC": "ovary",
                    "HCI-BCM": "breast",
                    "IRCC-CRC": "coadread", "IRCC-GC": "stomach", "LIH": "brain", "MDAnderson-CCH": "os", "NKI": "brca",
                    "SANG": "bowel",
                    "UMCG": "ovary", "UOC-BC": "brca", "UOM-BC": "brca", "VHIO-BC": "brca", "VHIO-CRC": "bowel",
                    "VHIO-PC": "pancreas"}
    if provider in cancer_types.keys():
        return cancer_types[provider]
    else:
        return "mixed"


def generate_meta_study_file(in_path, study, provider, out_path):
    df = pd.read_json("https://www.cancermodels.org/api/provider_group?abbreviation=eq." + provider)
    if provider != "DFCI-CPDM":
        provider_name_full = str(df["name"][0])
    elif provider == "DFCI-CPDM":
        provider_name_full = "Center for Patient Derived Models, Dana-Farber Cancer Institute"
    cancer_type = get_study_cancer_type(provider)
    meta_study_df = pd.DataFrame(columns=[0, 1])
    meta_study_df.loc[0, :] = ["type_of_cancer:", cancer_type]  # "mixed"]
    meta_study_df.loc[1, :] = ["cancer_study_identifier:", study]
    meta_study_df.loc[2, :] = ["name:", provider_name_full]
    meta_study_df.loc[3, :] = ["description:", get_provider_description(provider)]
    meta_study_df.loc[4, :] = ["groups:", "PUBLIC"]
    # meta_study_df.loc[5,:] = ["short_name:", "PDCMs ("+provider+")"]
    meta_study_df.loc[5, :] = ["reference_genome:", "hg38"]
    # meta_study_df.loc[6,:] = ["add_global_case_list:", "true"]
    meta_study_df.to_csv(join(out_path, "meta_study.txt"), sep="\t", index=False, header=False)


def generate_meta_clinical_files(study, out_path, type):
    meta_df = pd.DataFrame(columns=[0, 1])
    meta_df.loc[0, :] = ["cancer_study_identifier:", study]
    meta_df.loc[1, :] = ["genetic_alteration_type:", "CLINICAL"]
    meta_df.loc[2, :] = ["datatype:", type.upper() + "_ATTRIBUTES"]
    meta_df.loc[3, :] = ["data_filename:", "data_clinical_" + type + ".txt"]
    meta_df.to_csv(join(out_path, "meta_clinical_" + type + ".txt"), sep="\t", index=False, header=False)


def generate_meta_files(in_path, out_path, provider, study):
    generate_meta_study_file(in_path, study, provider, out_path)
    generate_meta_clinical_files(study, out_path, "sample")
    generate_meta_clinical_files(study, out_path, "patient")


def generate_case_lists_script(conf, out_path, study):
    if not exists(join(out_path, "case_lists")):
        makedirs(join(out_path, "case_lists"))
    # command = "python3 {} -c {} -d {} -s {} -i {}".format(script, conf, join(out_path, "case_lists"), out_path, study)
    main(conf, join(out_path, "case_lists"), out_path, study, False, False)


def generate_c_bio_portal_files(in_path, out_path, provider):
    # case_script = "/Users/tushar/CancerModels/utils/cbioportal/datahub-study-curation-tools/generate-case-lists/generate_case_lists.py"
    case_conf = "/Users/tushar/CancerModels/utils/cbioportal/datahub-study-curation-tools/generate-case-lists/case_list_conf.txt"
    study = provider
    out_path = join(out_path, study)
    if not exists(out_path):
        makedirs(out_path)
    generate_meta_files(in_path, out_path, provider, study)
    generate_clinical_sample(in_path, out_path, provider)
    generate_clinical_patient(in_path, out_path, provider)
    if provider != "CRL" and not provider.__contains__("Curie"):
        generate_timeline_data(in_path, out_path, study, provider)
        if exists(join(in_path, "mut")):
            generate_mutation_files(in_path, study, provider, out_path)
        if exists(join(in_path, "expression")):
            generate_expression_files(in_path, study, provider, out_path)
        if exists(join(in_path, "cna")):
            generate_cna_files(in_path, study, provider, out_path)
    generate_case_lists_script(case_conf, out_path, study)


# %% md
### Common env
# %%
start_dir = getcwd()
reference_df = get_hugo2ncbi()

if len(sys.argv)>0:
    home = sys.argv[1]
    out_path = sys.argv[2]
    providers = sorted(get_dirs(home))

    for i in tqdm(range(0, len(providers)),
                  desc="Generating cBioPortal data: "):  ## get_dirs will get the provider dirs in updog
        provider = providers[i]
        generate_c_bio_portal_files(join(home, provider), out_path, provider)
