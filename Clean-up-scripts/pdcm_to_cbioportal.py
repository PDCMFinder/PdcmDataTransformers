#!/usr/bin/env python
import sys
from os import listdir, getcwd, makedirs, remove
from os.path import isfile, join, isdir, exists
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import cpu_count


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


def get_mappings(home):
    mappings_path = home.replace('data/UPDOG/', '')
    mappings_path = join(mappings_path, 'mapping/diagnosis_mappings.json')
    mappings = pd.read_json(mappings_path)[['mappingValues', 'mappedTermLabel']]
    mappings = pd.concat([mappings, pd.json_normalize(mappings['mappingValues'])], axis=1)
    mappings = mappings.drop('mappingValues', axis=1).sort_values(by='DataSource').reset_index(drop=True)
    df = pd.read_json(
        'https://dev.cancermodels.org/api/model_metadata?select=cancer_system,histology').drop_duplicates().reset_index(
        drop=True)
    mappings = mappings.merge(df, left_on='mappedTermLabel', right_on='histology', how='left').fillna('-')
    return mappings


class cbioportal_case_lists:
    def __init__(self, filename, study_dir, study, overwrite, verbose):
        self.config_file_name = filename
        self.case_list_dir = join(study_dir, "case_lists")
        if not exists(self.case_list_dir):
            makedirs(self.case_list_dir)
        self.study_dir = study_dir
        self.study_id = study
        self.overwrite = overwrite
        self.verbose = verbose
        self.CASE_LIST_CONFIG_HEADER_COLUMNS = ["CASE_LIST_FILENAME", "STAGING_FILENAME", "META_STABLE_ID",
                                                "META_CASE_LIST_CATEGORY", "META_CANCER_STUDY_ID",
                                                "META_CASE_LIST_NAME",
                                                "META_CASE_LIST_DESCRIPTION"]
        self.CASE_LIST_UNION_DELIMITER = "|"
        self.CASE_LIST_INTERSECTION_DELIMITER = "&"
        self.MUTATION_STAGING_GENERAL_PREFIX = "data_mutations"
        self.SEQUENCED_SAMPLES_FILENAME = "sequenced_samples.txt"
        self.MUTATION_CASE_LIST_META_HEADER = "sequenced_samples"
        self.MUTATION_CASE_ID_COLUMN_HEADER = "Tumor_Sample_Barcode"
        self.SAMPLE_ID_COLUMN_HEADER = "SAMPLE_ID"
        self.NON_CASE_IDS = frozenset(
            ["MIRNA", "LOCUS", "ID", "GENE SYMBOL", "ENTREZ_GENE_ID", "HUGO_SYMBOL", "LOCUS ID", "CYTOBAND",
             "COMPOSITE.ELEMENT.REF", "HYBRIDIZATION REF"])
        self.CANCER_STUDY_TAG = "<CANCER_STUDY>"
        self.NUM_CASES_TAG = "<NUM_CASES>"

    def generate_case_lists(self):
        header = []
        with open(self.config_file_name, 'r') as case_list_config_file:
            # get header and validate
            header = case_list_config_file.readline().rstrip('\n').rstrip('\r').split('\t')
            # check full header matches what we expect
            for column in self.CASE_LIST_CONFIG_HEADER_COLUMNS:
                if column not in header:
                    print >> sys.stderr, "ERROR: column '%s' is not in '%s'" % (column, self.config_file_name)
                    sys.exit(2)

            for line in case_list_config_file:
                line = line.rstrip('\n').rstrip('\r')
                config_fields = line.split('\t')
                case_list_filename = config_fields[header.index("CASE_LIST_FILENAME")]
                staging_filename_list = config_fields[header.index("STAGING_FILENAME")]
                case_list_file_full_path = join(self.case_list_dir, case_list_filename)
                if isfile(case_list_file_full_path) and not self.overwrite:
                    if self.verbose:
                        print("LOG: generate_case_lists(), '%s' exists and overwrite is false, skipping caselist..." % (
                            case_list_filename))
                    continue

                # might be single staging file
                staging_filenames = []
                # union (like all cases)
                union_case_list = self.CASE_LIST_UNION_DELIMITER in staging_filename_list
                # intersection (like complete or cna-seq)
                intersection_case_list = self.CASE_LIST_INTERSECTION_DELIMITER in staging_filename_list
                delimiter = self.CASE_LIST_UNION_DELIMITER if union_case_list else self.CASE_LIST_INTERSECTION_DELIMITER
                staging_filenames = staging_filename_list.split(delimiter)
                if self.verbose:
                    print("LOG: generate_case_lists(), staging filenames: %s" % (",".join(staging_filenames)))

                # if this is intersection all staging files must exist
                if intersection_case_list and \
                        not all([isfile(join(self.study_dir, intersection_filename)) for intersection_filename in
                                 staging_filenames]):
                    continue

                # this is the set we will pass to write_case_list_file
                case_set = set([])
                # this indicates the number of staging files processed -
                # used to verify that an intersection should be written
                num_staging_files_processed = 0
                for staging_filename in staging_filenames:
                    if self.verbose:
                        print("LOG: generate_case_lists(), processing staging file '%s'" % (staging_filename))
                    # compute the case set
                    case_list = []
                    case_list = self.get_case_list_from_staging_file(staging_filename)

                    if len(case_list) == 0:
                        if self.verbose:
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
                    if self.verbose:
                        print("LOG: generate_case_lists(), calling write_case_list_file()...")

                    # do not write out complete cases file unless we've processed all the files required
                    if intersection_case_list and num_staging_files_processed != len(staging_filenames):
                        if self.verbose:
                            print(
                                "LOG: generate_case_lists(), number of staging files processed (%d) != number of staging files required (%d) for '%s', skipping call to write_case_list_file()..." % (
                                    num_staging_files_processed, len(staging_filenames), case_list_filename))
                    else:
                        self.write_case_list_file(header, config_fields, case_list_file_full_path, case_set)
                elif self.verbose:
                    print(
                        "LOG: generate_case_lists(), case_set.size() == 0, skipping call to write_case_list_file()...")

    def get_case_list_from_staging_file(self, staging_filename):
        if self.verbose:
            print("LOG: get_case_list_from_staging_file(), '%s'" % (staging_filename))

        case_set = set([])

        # if we are processing mutations data and a SEQUENCED_SAMPLES_FILENAME exists, use it
        if self.MUTATION_STAGING_GENERAL_PREFIX in staging_filename:
            sequenced_samples_full_path = join(self.study_dir, self.SEQUENCED_SAMPLES_FILENAME)
            if isfile(sequenced_samples_full_path):
                if self.verbose:
                    print(
                        "LOG: get_case_list_from_staging_file(), '%s' exists, calling get_case_list_from_sequenced_samples_file()" % (
                            self.SEQUENCED_SAMPLES_FILENAME))
                return self.get_case_list_from_sequenced_samples_file(sequenced_samples_full_path)

        staging_file_full_path = join(self.study_dir, staging_filename)
        if not isfile(staging_file_full_path):
            return []

        # staging file
        with open(staging_file_full_path, 'r') as staging_file:
            id_column_index = 0
            process_header = True
            for line in staging_file:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    if line.startswith('#' + self.MUTATION_CASE_LIST_META_HEADER + ':'):
                        # split will split on any whitespace, tabs or any number of consecutive spaces
                        return line[len(self.MUTATION_CASE_LIST_META_HEADER) + 2:].strip().split()
                    continue  # this is a comment line, skip it
                values = line.split('\t')

                # is this the header line?
                if process_header:
                    # look for MAF file case id column header
                    # if this is not a MAF file and header contains the case ids, return here
                    # we are assuming the header contains the case ids because SAMPLE_ID_COLUMN_HEADER is missing
                    if self.MUTATION_CASE_ID_COLUMN_HEADER not in values and self.SAMPLE_ID_COLUMN_HEADER not in [
                        x.upper() for x
                        in
                        values]:
                        if self.verbose:
                            print(
                                "LOG: get_case_list_from_staging_file(), this is not a MAF header but has no '%s' column, we assume it contains sample ids..." % (
                                    self.SAMPLE_ID_COLUMN_HEADER))
                        for potential_case_id in values:
                            # check to filter out column headers other than sample ids
                            if potential_case_id.upper() in self.NON_CASE_IDS:
                                continue
                            case_set.add(potential_case_id)
                        break  # got case ids from header, don't read the rest of the file
                    else:
                        # we know at this point one of these columns exists, so no fear of ValueError from index method
                        id_column_index = values.index(
                            self.MUTATION_CASE_ID_COLUMN_HEADER) if self.MUTATION_CASE_ID_COLUMN_HEADER in values else [
                            x.upper()
                            for
                            x in
                            values].index(
                            self.SAMPLE_ID_COLUMN_HEADER)
                        if self.verbose:
                            print(
                                "LOG: get_case_list_from_staging_file(), this is a MAF or clinical file, samples ids in column with index: %d" % (
                                    id_column_index))
                    process_header = False
                    continue  # done with header, move on to next line
                case_set.add(values[id_column_index])

        return list(case_set)

    def get_case_list_from_sequenced_samples_file(self, sequenced_samples_full_path):
        if self.verbose:
            print("LOG: get_case_list_from_sequenced_samples_file, '%s'", sequenced_samples_full_path)

        case_set = set([])
        with open(sequenced_samples_full_path, 'r') as sequenced_samples_file:
            for line in sequenced_samples_file:
                case_set.add(line.rstrip('\n'))

        if self.verbose:
            print("LOG: get_case_list_from_sequenced_samples_file, case set size: %d" % (len(case_set)))

        return list(case_set)

    def write_case_list_file(self, case_list_config_header, case_list_config_fields, case_list_full_path, case_set):
        if self.verbose:
            print("LOG: write_case_list_file(), '%s'" % (case_list_full_path))
        with open(case_list_full_path, 'w') as case_list_file:
            case_list_file.write("cancer_study_identifier: " + self.study_id + "\n")
            stable_id = case_list_config_fields[case_list_config_header.index("META_STABLE_ID")].replace(
                self.CANCER_STUDY_TAG,
                self.study_id)
            case_list_file.write("stable_id: " + stable_id + "\n")
            case_list_file.write(
                "case_list_name: " + case_list_config_fields[
                    case_list_config_header.index("META_CASE_LIST_NAME")] + "\n")
            case_list_description = case_list_config_fields[
                case_list_config_header.index("META_CASE_LIST_DESCRIPTION")].replace(self.NUM_CASES_TAG,
                                                                                     str(len(case_set)))
            case_list_file.write("case_list_description: " + case_list_description + "\n")
            case_list_file.write("case_list_category: " + case_list_config_fields[
                case_list_config_header.index("META_CASE_LIST_CATEGORY")] + "\n")
            case_list_file.write("case_list_ids: " + '\t'.join(case_set) + "\n")

    def main(self):
        if self.verbose:
            print("LOG: case_list_config_file='%s'" % (self.config_file_name))
            print("LOG: case_list_dir='%s'" % (self.case_list_dir))
            print("LOG: study_dir='%s'" % (self.study_dir))
            print("LOG: study_id='%s'" % (self.study_id))
            print("LOG: overwrite='%s'" % (self.overwrite))
            print("LOG: verbose='%s'" % (self.verbose))

        if not isfile(self.config_file_name):
            print("ERROR: case list configuration file '%s' does not exist or is not a file" % (
                self.config_file_name),
                  file=sys.stderr)
            sys.exit(2)

        if not isdir(self.case_list_dir):
            print("ERROR: case list file directory '%s' does not exist or is not a directory" % (self.case_list_dir),
                  file=sys.stderr)
            sys.exit(2)

        if not isdir(self.study_dir):
            print("ERROR: study directory '%s' does not exist or is not a directory" % (self.study_dir),
                  file=sys.stderr)
            sys.exit(2)

        self.generate_case_lists()


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
                    "IRCC-CRC": "coadread", "IRCC-GC": "stomach", "LIH": "brain", "MDAnderson-CCH": "os",
                    "NKI": "brca",
                    "SANG": "bowel",
                    "UMCG": "ovary", "UOC-BC": "brca", "UOM-BC": "brca", "VHIO-BC": "brca", "VHIO-CRC": "bowel",
                    "VHIO-PC": "pancreas"}
    if provider in cancer_types.keys():
        return cancer_types[provider]
    else:
        return "mixed"


def write_clinical_file(df, path, file_type):
    df.to_csv(join(path, "data_clinical_" + file_type + ".txt"), sep="\t", index=False)


def add_headers_clinical_patient():
    # Column headers - The attribute Display Names: The display name for each clinical attribute
    column_headers = ["#Patient Identifier", "Sex", "Diagnosis Age", "Overall Survival (Months)",
                      "Overall Survival Status"]
    # Column descriptions - The attribute Descriptions: Long(er) description of each clinical attribute
    column_description = ["#Identifier to uniquely specify a patient.", "Sex",
                          "Age at which a condition or disease was first diagnosed.",
                          "Overall survival in months since initial diagonosis.",
                          "Overall patient survival status."]
    # Column data type - The attribute Datatype: The datatype of each clinical attribute (must be one of: STRING, NUMBER, BOOLEAN)
    column_data_type = ["#STRING", "STRING", "NUMBER", "NUMBER", "STRING"]
    # Column priority - A number which indicates the importance of each attribute. A higher number indicates a higher priority.
    column_priority = ["#1", "1", "1", "1", "1"]
    # Column headers for validation
    column_header_validation = ["PATIENT_ID", "SEX", "AGE", "OS_MONTHS", "OS_STATUS"]
    return [column_headers, column_description, column_data_type, column_priority, column_header_validation]


def add_headers_clinical_sample():
    # Column headers - The attribute Display Names: The display name for each clinical attribute
    column_headers = ["#Patient Identifier", "Sample Identifier", 'Sample Origin', 'Passage', "Tumor Type",
                      "Cancer Type", "Cancer Type Detailed",
                      "Primary site", "Tumor Grade", "Model type", "Model ID"]
    # Column descriptions - The attribute Descriptions: Long(er) description of each clinical attribute
    column_description = ["#Identifier to uniquely specify a patient.", "A unique sample identifier.",
                          "Sample Origin", 'Passage',
                          "The type of tumour sample (i.e., normal, primary, met, recurrence).", "Cancer Type",
                          "Cancer Type Detailed",
                          "Site of the primary tumor where primary cancer is originating from (may not correspond to the site of the current tissue sample). No abbreviations.",
                          "The implanded tumour grade value.", "Type of patient derived cancer model",
                          "Unique identifier for the PDCMs"]
    # Column data type - The attribute Datatype: The datatype of each clinical attribute (must be one of: STRING, NUMBER, BOOLEAN)
    column_data_type = ["#STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING",
                        "STRING", "STRING"]
    # Column priority - A number which indicates the importance of each attribute. A higher number indicates a higher priority.
    column_priority = ["#1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"]
    # Column headers for validation
    column_header_validation = ["PATIENT_ID", "SAMPLE_ID", "SAMPLE_ORIGIN", "PASSAGE", "SAMPLE_TYPE", "CANCER_TYPE",
                                "CANCER_TYPE_DETAILED",
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


def get_model_type(df, path, provider):
    pdx, cell = pd.DataFrame(), pd.DataFrame()
    if exists(join(path, provider + '_metadata-pdx_model.tsv')):
        pdx = read_metadata_without_fields(join(path, provider + '_metadata-pdx_model.tsv'))
        pdx['model_type'] = 'PDX'
        pdx = pdx[['model_id', 'model_type']]

    if exists(join(path, provider + '_metadata-cell_model.tsv')):
        cell = read_metadata_without_fields(join(path, provider + '_metadata-cell_model.tsv'))[['model_id', 'type']]
        cell['model_type'] = cell['type']
        cell = cell[['model_id', 'model_type']]
    if pdx.shape[0] > 0 and cell.shape[0] > 0:
        temp = pd.concat([pdx, cell]).reset_index(drop=True)
    elif pdx.shape[0] > 0 and cell.shape[0] == 0:
        temp = pdx
    elif pdx.shape[0] == 0 and cell.shape[0] > 0:
        temp = cell
    return df.merge(temp, on='model_id', how='left')


def get_samples_mms(df, col, path, provider):
    df['Sample Origin'] = ''
    df['Passage'] = ''
    if not exists(join(path, provider + '_molecular_metadata-sample.tsv')):
        return df
    mms = read_metadata_without_fields(join(path, provider + '_molecular_metadata-sample.tsv'))
    missing_samples = list(set(mms['sample_id'].unique()) - set(df['Sample Identifier'].unique()))
    missing = mms[mms['sample_id'].isin(missing_samples)][['model_id', 'sample_id', 'sample_origin', 'passage']]
    missing = missing[missing['sample_origin'].str.lower() != 'patient']
    missing['Sample Origin'] = missing['sample_origin']
    missing['Passage'] = missing['passage']
    temp = missing.merge(df, right_on='Model ID', left_on='model_id', how='left')
    temp['Sample Identifier'] = temp['sample_id']
    temp['Sample Origin'] = temp['Sample Origin_x']
    temp['Passage'] = temp['Passage_x']
    temp = temp[col]
    return pd.concat([df, temp]).reset_index(drop=True)[col]


def filter_sample_ids(x, sample_ids):
    out = [val for val in sample_ids if val in x]
    if len(out) > 0 and out[0] != '':
        return out[0]
    else:
        "No match"


def filter_samples(df, col, out_path):
    sample_ids = list(
        pd.read_csv(join(out_path, "data_clinical_sample.txt"), sep="\t")["Sample Identifier"].astype(str))[4:]
    df[col] = df[col].astype(str)
    df = df[df[col].apply(lambda x: any(str(val) in x for val in sample_ids))]
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


class cBioPortal:
    def __init__(self, home, output_dir):
        self.study, self.provider_path, self.study_dir, self.mms = None, None, None, None
        self.home = home
        self.providers = sorted(get_dirs(home))
        self.mappings = get_mappings(home)
        self.output_dir = output_dir
        self.case_conf = "case_list_conf.txt"
        self.threads = int(cpu_count() / 2)

    def get_platform(self, type):
        df = self.mms[self.mms["molecular_characterisation_type"] == type]
        df["platform_name"] = df["library_strategy"] + " -  " + df["instrument_model"]
        if self.study == "JAX":
            if type == "mutation":
                return "WES"
            if type == "expression":
                return "RNA-Seq"
            if type == "copy number alteration":
                return "SNP"
        return ", ".join(df["platform_name"])

    def generate_meta_cna(self, platform, datatype):
        meta_df = pd.DataFrame(columns=[0, 1])
        meta_df.loc[0, :] = ["cancer_study_identifier:", self.study]
        meta_df.loc[1, :] = ["genetic_alteration_type:", "COPY_NUMBER_ALTERATION"]
        if datatype == "log":
            meta_df.loc[2, :] = ["datatype:", "LOG2-VALUE"]
            meta_df.loc[3, :] = ["stable_id:", "log2CNA"]
            out_name = "log2_cna"
            meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "false"]
            meta_df.loc[5, :] = ["profile_description:", "Log2 copy-number data from " + platform]
            meta_df.loc[6, :] = ["profile_name:", "Log2 copy-number values"]
            meta_df.loc[7, :] = ["data_filename:", "data_log2_cna.txt"]

        elif datatype == "gistic":
            meta_df.loc[2, :] = ["datatype:", "DISCRETE"]
            meta_df.loc[3, :] = ["stable_id:", "gistic"]
            out_name = "cna"
            meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "false"]
            meta_df.loc[5, :] = ["profile_description:",
                                 "Putative copy-number from GISTIC 2.0. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification. " + platform]
            meta_df.loc[6, :] = ["profile_name:", "Putative copy-number alterations from GISTIC"]
            meta_df.loc[7, :] = ["data_filename:", "data_cna.txt"]
        meta_df.to_csv(join(self.study_dir, "meta_" + out_name + ".txt"), sep="\t", index=False, header=False)

    def generate_cna_data(self, platform):
        df = read_mol_data(join(self.provider_path, "cna"))[["sample_id", "symbol", "log2r_cna", "gistic_value"]]
        col_value = "log2r_cna"
        datatype = "log"
        out_name = "log2_cna"
        df = filter_samples(df, "sample_id", self.study_dir)
        if len(df) == 0:
            return None
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
            self.generate_meta_cna(platform, datatype)
            df[out_cols].to_csv(join(self.study_dir, "data_" + out_name + ".txt"), index=False, sep="\t")

    def generate_meta_expression(self, platform, datatype):
        meta_df = pd.DataFrame(columns=[0, 1])
        meta_df.loc[0, :] = ["cancer_study_identifier:", self.study]
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
        meta_df.to_csv(join(self.study_dir, "meta_mrna_seq_rpkm.txt"), sep="\t", index=False, header=False)

    def generate_expression_data(self, platform):
        df = read_mol_data(join(self.provider_path, "expression"))[["sample_id", "symbol", "rnaseq_fpkm", "z_score"]]
        df = filter_samples(df, "sample_id", self.study_dir)
        if len(df) == 0:
            return None
        df["Hugo_Symbol"] = df["symbol"]
        no_fpkm = df['rnaseq_fpkm'].isna().all()
        if no_fpkm and df['z_score'].isna().all():
            return None
        datatype = "mrna"
        value_column = "rnaseq_fpkm"
        if no_fpkm:
            print("Using Z score: " + self.study)
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
            self.generate_meta_expression(platform, datatype)
            df[out_cols].to_csv(join(self.study_dir, "data_mrna_seq_rpkm.txt"), sep="\t", index=False)

    def generate_meta_mutation(self, platform):
        meta_df = pd.DataFrame(columns=[0, 1])
        meta_df.loc[0, :] = ["cancer_study_identifier:", self.study]
        meta_df.loc[1, :] = ["genetic_alteration_type:", "MUTATION_EXTENDED"]
        meta_df.loc[2, :] = ["datatype:", "MAF"]
        meta_df.loc[3, :] = ["stable_id:", "mutations"]
        meta_df.loc[4, :] = ["show_profile_in_analysis_tab:", "true"]
        meta_df.loc[5, :] = ["profile_description:", "Mutation data from " + platform]
        meta_df.loc[6, :] = ["profile_name:", "Mutations"]
        meta_df.loc[7, :] = ["variant_classification_filter:",
                             "De_novo_Start_InFrame, De_novo_Start_OutOfFrame, Unknown"]
        meta_df.loc[8, :] = ["data_filename:", "data_mutations_extended.txt"]
        meta_df.to_csv(join(self.study_dir, "meta_mutations_extended.txt"), sep="\t", index=False, header=False)

    def generate_mutation_data(self, platform):
        mapper = {"sample_id": "Tumor_Sample_Barcode", "seq_start_position": "Start_Position",
                  "chromosome": "Chromosome",
                  "ref_allele": "Reference_Allele", "alt_allele": "Tumor_Seq_Allele2", "symbol": "Hugo_Symbol",
                  "consequence": "Consequence", "amino_acid_change": "HGVSp_Short",
                  "ensembl_gene_id": "Gene", "ensembl_transcript_id": "Transcript_ID", "codon_change": "Codons",
                  "ncbi_transcript_id": "RefSeq"}
        mut_df = read_mol_data(join(self.provider_path, "mut")).rename(columns=mapper).fillna("")
        mut_df = filter_samples(mut_df, "Tumor_Sample_Barcode", self.study_dir)
        if len(mut_df) == 0:
            return None
        out_column = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position",
                      "End_Position", "Strand", "Consequence", "Variant_Classification", "Variant_Type",
                      "Reference_Allele",
                      "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode",
                      "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
                      "Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
                      "Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status",
                      "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", "BAM_File", "Sequencer",
                      "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count", "HGVSc", "HGVSp", "HGVSp_Short",
                      "Transcript_ID", "RefSeq", "Protein_position", "Codons", "Hotspot"]
        mut_df["End_Position"] = pd.to_numeric(mut_df["Start_Position"]) + np.where(
            mut_df.Reference_Allele.str.len() >= mut_df.Tumor_Seq_Allele2.str.len(), mut_df.Reference_Allele.str.len(),
            mut_df.Tumor_Seq_Allele2.str.len()) - 1
        mut_df["NCBI_Build"] = "GRCh38"
        mut_df["Strand"] = "+"
        mut_df["HGVSc"] = mut_df["Transcript_ID"] + "c." + mut_df["coding_sequence_change"]
        mut_df['Protein_position'] = np.where(mut_df['HGVSp_Short'].fillna('') != '', mut_df['HGVSp_Short'].str[1:-1],
                                              '')
        mut_df["HGVSp_Short"] = np.where(mut_df["HGVSp_Short"].fillna('') != "",
                                         "p." + mut_df["HGVSp_Short"].str.replace('.0', ''), "")
        mut_df["Entrez_Gene_Id"] = mut_df['Hugo_Symbol'].apply(map_gene_symbol_to_id)
        mut_df = mut_df[mut_df["Entrez_Gene_Id"] != ""]
        mut_df['Hugo_Symbol'] = mut_df[["Entrez_Gene_Id", 'Hugo_Symbol']].apply(lambda x: map_id_to_hugo(x[0], x[1]),
                                                                                axis=1)
        mut_df['RefSeq'] = [
            "" if egid == "" else reference_df.loc[reference_df['NCBI Gene ID'] == egid, 'RefSeq IDs'].values[0] for
            egid in mut_df['Entrez_Gene_Id']]
        # mut_df["User_Amino_Acid_Change"] = mut_df["HGVSp_Short"]

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
        mut_df = mut_df.apply(lambda x: handle_frameshift(x), axis=1)
        mut_df["Variant_Classification"] = mut_df["Consequence"].apply(lambda x: mapper.get(x, "Unknown"))
        mut_df["Variant_Classification"] = mut_df["Variant_Classification"].replace(mapper)
        mut_df["Variant_Type"] = mut_df["variant_class"].replace("SNV", "SNP").replace("insertion", "INS").replace(
            "deletion", "DEL")
        for column in out_column:
            if column not in mut_df.columns:
                mut_df[column] = ""
        if mut_df.shape[0] > 0:
            self.generate_meta_mutation(platform)
            mut_df[out_column].to_csv(join(self.study_dir, "data_mutations_extended.txt"), index=False, sep="\t")

    def generate_molecular_data_files(self, dt):
        self.mms = read_metadata_without_fields(
            join(self.provider_path, self.study + "_molecular_metadata-platform.tsv")).fillna("")
        platform = self.get_platform(dt)
        if dt == 'mut' and exists(join(self.provider_path, "mut")):
            self.generate_mutation_data(platform)
        elif dt == 'expression' and exists(join(self.provider_path, "expression")):
            self.generate_expression_data(platform)
        elif dt == 'copy number alteration' and exists(join(self.provider_path, "cna")):
            self.generate_cna_data(platform)
        elif dt == 'timeline':
            self.generate_timeline_data()

    def generate_meta_timeline(self, type):
        meta_df = pd.DataFrame(columns=[0, 1])
        meta_df.loc[0, :] = ["cancer_study_identifier:", self.study]
        meta_df.loc[1, :] = ["genetic_alteration_type:", "CLINICAL"]
        meta_df.loc[2, :] = ["datatype:", "TIMELINE"]
        meta_df.loc[7, :] = ["data_filename:", "data_timeline_" + type + ".txt"]
        meta_df.to_csv(join(self.study_dir, "meta_timeline_" + type + ".txt"), sep="\t", index=False, header=False)

    def generate_timeline_data_file(self, type):
        ps = read_metadata_without_fields(join(self.provider_path, self.study + "_metadata-patient_sample.tsv"))
        self.generate_meta_timeline(type)
        if type == 'specimen':
            out_cols = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SAMPLE_ID", "SPECIMEN_SITE", "SOURCE"]
            ps["PATIENT_ID"], ps["SAMPLE_ID"], ps["STOP_DATE"] = ps["patient_id"], ps["sample_id"], ""
            ps['START_DATE'], ps['SOURCE'], ps['EVENT_TYPE'] = 0, self.study, type.upper()
            ps["SPECIMEN_SITE"] = ps["collection_site"]
            ids = list(pd.read_csv(join(self.study_dir, "data_clinical_patient.txt"), sep="\t")["#Patient Identifier"])
            ps = ps[ps["patient_id"].isin(ids)]
            ps[out_cols].to_csv(join(self.study_dir, "data_timeline_" + type + ".txt"), sep="\t", index=False)
        if type == "treatment":
            treatment = read_metadata_with_fields(
                join(self.provider_path, "treatment", self.study + "_patienttreatment-Sheet1.tsv"))
            dates = ps[["patient_id", "collection_date"]]
            try:
                dates["date_time"] = pd.to_datetime(
                    dates["collection_date"].replace("(?i)not provided", "", regex=True).replace("(?i)not collected",
                                                                                                 "",
                                                                                                 regex=True).str.replace(
                        " ", "-"), format='%b-%y')
            except:
                try:
                    dates["date_time"] = pd.to_datetime(
                        dates["collection_date"].replace("(?i)not provided", "", regex=True).replace(
                            "(?i)not collected",
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
                    treatment["treatment_starting_date"].str.replace('Not provided', '',
                                                                     flags=re.IGNORECASE).str.replace(
                        'Not collected', '', flags=re.IGNORECASE).str.replace(" ", "-"), format='%b-%y')
            except:
                try:
                    treatment["start_date_time"] = pd.to_datetime(
                        treatment["treatment_starting_date"].str.replace('Not provided', '',
                                                                         flags=re.IGNORECASE).str.replace(
                            'Not collected',
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
            try:
                treatment["START_DATE"] = (treatment['start_date_time'] - treatment['date_time']).dt.days
            except:
                treatment["START_DATE"] = 0
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

            treatment[out_cols].to_csv(join(self.study_dir, "data_timeline_" + type + ".txt"), sep="\t", index=False)
        if type == "lab_test":
            cyto = read_metadata_with_fields(
                join(self.provider_path, "biomarker", self.study + "_biomarker-Sheet1.tsv"))
            cyto = cyto.merge(ps[["patient_id", "sample_id"]], on="sample_id", how="inner")
            out_cols = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT"]
            cyto["PATIENT_ID"], cyto["START_DATE"], cyto["STOP_DATE"] = cyto["patient_id"], 0, ""
            cyto["EVENT_TYPE"] = type.upper()
            cyto["TEST"], cyto["RESULT"] = cyto['symbol'], cyto['marker_status']
            cyto[out_cols].to_csv(join(self.study_dir, "data_timeline_" + type + ".txt"), sep="\t", index=False)

    def generate_timeline_data(self):
        self.generate_timeline_data_file('specimen')
        if exists(join(self.provider_path, "treatment")):
            self.generate_timeline_data_file('treatment')
        if exists(join(self.provider_path, "biomarker")):
            self.generate_timeline_data_file('lab_test')

    def clinical_patient_data_transformation(self):
        out_df = generate_clinical_df(add_headers_clinical_patient())
        mapper = {"patient_id": "#Patient Identifier", "sex": "Sex", "age_at_initial_diagnosis": "Diagnosis Age"}
        temp = read_metadata_without_fields(join(self.provider_path, self.study + "_metadata-patient.tsv")).rename(
            columns=mapper)

        temp["Overall Survival (Months)"] = ""
        temp["Overall Survival Status"] = ""
        cs = pd.read_csv(join(self.study_dir, "data_clinical_sample.txt"), sep="\t")
        patient_ids = list(cs["#Patient Identifier"])
        temp = temp[temp["#Patient Identifier"].isin(patient_ids)]
        out_df = pd.concat([out_df, temp[out_df.columns]], ignore_index=True)

        return out_df.replace("(?i)not provided", "", regex=True).replace("(?i)not collected", "",
                                                                          regex=True).drop_duplicates(
            subset=["#Patient Identifier"])

    def clinical_patient_sample_data_transformation(self):
        out_df = generate_clinical_df(add_headers_clinical_sample())
        mapper = {"patient_id": "#Patient Identifier", "sample_id": "Sample Identifier", "tumour_type": "Tumor Type OG",
                  "primary_site": "Primary site OG", "grade": "Tumor Grade"}
        temp = read_metadata_without_fields(
            join(self.provider_path, self.study + "_metadata-patient_sample.tsv")).rename(columns=mapper)
        temp = get_model_type(temp, self.provider_path, self.study)
        temp['diagnosis'] = temp['diagnosis'].str.lower()
        temp['tumortype'] = temp['Tumor Type OG'].str.lower()
        temp['origin'] = temp['Primary site OG'].str.lower()
        temp = temp.merge(self.mappings[self.mappings['DataSource'] == self.study.lower()],
                          left_on=['diagnosis', 'tumortype', 'origin'],
                          right_on=['SampleDiagnosis', 'TumorType', 'OriginTissue'], how='left')
        # temp = temp.merge(get_meta_from_api(provider), on="model_id", how="inner").reset_index(drop=True)
        temp["Cancer Type"] = temp["cancer_system"]
        temp["Cancer Type Detailed"] = temp["mappedTermLabel"].str.title()
        temp["Model type"] = temp["model_type"]
        temp["Tumor Type"] = temp["tumortype"].str.title()
        temp["Primary site"] = temp["origin"].str.title()
        temp["Model ID"] = temp["model_id"]
        temp["Tumor Grade"] = temp["Tumor Grade"].str.replace('Not provided', '', flags=re.IGNORECASE).str.replace(
            'Not collected', '', flags=re.IGNORECASE).replace("p.*", "", regex=True).replace("T.*", "",
                                                                                             regex=True).replace(
            " \(.*", "", regex=True).replace(";.*", "", regex=True).str.replace('Not collected', '',
                                                                                flags=re.IGNORECASE).replace("Moderate",
                                                                                                             "moderately differentiated",
                                                                                                             regex=True).replace(
            "Poorly", "poorly differentiated", regex=True).replace("Well", "well differentiated", regex=True).replace(
            "High-Grade", "high", regex=True).replace("High grade", "high", regex=True).replace("Grade ", "G",
                                                                                                regex=True)
        temp = get_samples_mms(temp, out_df.columns, self.provider_path, self.study)
        out_df = pd.concat([out_df, temp], ignore_index=True)
        return out_df.replace("(?i)not provided", "", regex=True).replace("(?i)not collected", "",
                                                                          regex=True).drop_duplicates(
            subset=["#Patient Identifier", "Sample Identifier"])

    def generate_meta_study_file(self):
        df = pd.read_json("https://www.cancermodels.org/api/provider_group?abbreviation=eq." + self.study)
        if self.study != "DFCI-CPDM":
            provider_name_full = str(df["name"][0])
        elif self.study == "DFCI-CPDM":
            provider_name_full = "Center for Patient Derived Models, Dana-Farber Cancer Institute"
        cancer_type = get_study_cancer_type(self.study)
        meta_study_df = pd.DataFrame(columns=[0, 1])
        meta_study_df.loc[0, :] = ["type_of_cancer:", cancer_type]  # "mixed"]
        meta_study_df.loc[1, :] = ["cancer_study_identifier:", self.study]
        meta_study_df.loc[2, :] = ["name:", provider_name_full]
        meta_study_df.loc[3, :] = ["description:", get_provider_description(self.study)]
        meta_study_df.loc[4, :] = ["groups:", "PUBLIC"]
        # meta_study_df.loc[5,:] = ["short_name:", "PDCMs ("+provider+")"]
        meta_study_df.loc[5, :] = ["reference_genome:", "hg38"]
        # meta_study_df.loc[6,:] = ["add_global_case_list:", "true"]
        meta_study_df.to_csv(join(self.study_dir, "meta_study.txt"), sep="\t", index=False, header=False)

    def generate_meta_clinical_files(self, type):
        meta_df = pd.DataFrame(columns=[0, 1])
        meta_df.loc[0, :] = ["cancer_study_identifier:", self.study]
        meta_df.loc[1, :] = ["genetic_alteration_type:", "CLINICAL"]
        meta_df.loc[2, :] = ["datatype:", type.upper() + "_ATTRIBUTES"]
        meta_df.loc[3, :] = ["data_filename:", "data_clinical_" + type + ".txt"]
        meta_df.to_csv(join(self.study_dir, "meta_clinical_" + type + ".txt"), sep="\t", index=False, header=False)

    def meta_and_clinical_files(self):
        self.generate_meta_study_file()
        self.generate_meta_clinical_files('sample')
        self.generate_meta_clinical_files('patient')
        write_clinical_file(self.clinical_patient_sample_data_transformation(), self.study_dir, "sample")
        write_clinical_file(self.clinical_patient_data_transformation(), self.study_dir, "patient")

    def generate_c_bio_portal_files(self, provider):
        self.study = provider
        self.provider_path = join(self.home, provider)  # in
        self.study_dir = join(self.output_dir, provider)  # out
        dt = ['mut', 'copy number alteration', 'expression', 'timeline']
        if not exists(self.study_dir):
            makedirs(self.study_dir)
        self.meta_and_clinical_files()
        if provider != "CRL" and not provider.__contains__("Curie"):
            with ThreadPoolExecutor(max_workers=self.threads) as executor:
                executor.map(self.generate_molecular_data_files, dt)

        cases = cbioportal_case_lists(self.case_conf, self.study_dir, provider, False, False)
        cases.main()

    def main(self):
        for i in tqdm(range(0, len(self.providers)), desc="Generating cBioPortal data: "):
            self.generate_c_bio_portal_files(self.providers[i])


reference_df = get_hugo2ncbi()


def run(args):
    # Path to data # output path
    cBioPortal(args[1], args[2]).main()


if len(sys.argv) > 0:
    run(sys.argv)