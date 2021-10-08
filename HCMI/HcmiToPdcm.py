import pandas as pd

def align_and_append_sheet(template, other_sheet):
    alignment = template.align(other_sheet, join="outer")[1]
    return template.append(alignment).dropna(axis=0, how='all')

def generate_patient_id(model_name):
    patient_id = f'{model_name}_patient'
    suffix = model_name[-2:]
    if suffix == "-B" or suffix == "-A":
        patient_id = f'{model_name[:-2]}_patient'
    return patient_id

def generate_sample_id(model_name):
    return f'{model_name}_sample'

def cat_str_if_not_NA(preprend_str, str_to_check):
    built_str = ''
    if str_to_check != '--':
        built_str = preprend_str + str_to_check
    return built_str

def buildHistoryString(row):
    vital_status = row["Vital Status"]
    disease_status = row["Disease Status"]
    built_str = ''
    vital_status_str = cat_str_if_not_NA("vital_status:", vital_status)
    disease_status_str = cat_str_if_not_NA("disease_status:", disease_status)
    if vital_status_str != '':
        built_str = f'{vital_status_str},{disease_status_str}'
    return built_str

def generate_patient_sheet(hcmiDf, patient_template):
    patient_columns = ["Gender","Race", "Age At Diagnosis (Years)", "Name", "Vital Status", "Disease Status" ]
    column_name_map = {"Gender": "sex", "Race" : "ethnicity","Age At Diagnosis (Years)" :"age_at_initial_diagnosis", "Name": "patient_id"}
    HCMI_patient_columns = hcmiDf.loc[:, patient_columns].rename(columns=column_name_map)
    HCMI_patient_columns['history'] = hcmiDf.apply(buildHistoryString,axis=1)
    HCMI_patient_columns["patient_id"] = HCMI_patient_columns["patient_id"].apply(generate_patient_id)
    HCMI_patient_columns = HCMI_patient_columns.drop(columns=["Vital Status", "Disease Status"])
    align_and_append_sheet(patient_template, HCMI_patient_columns).to_csv("HCMI_metadata-patient.tsv", sep='\t', index=False)

def generate_sample_sheet(hcmiDf, sample_template, clinicalDf):
    sample_columns = ["Name", "Primary Site", "Clinical Tumor Diagnosis", "Histological Subtype", "Tissue Status", "Acquisition Site", "Age At Acquisition (Years)", "TNM Stage", "Histological Grade"]
    sample_name_map = { "Name": "model_id", "Primary Site" : "primary_site", "Tissue Status" : "tumour_type", "Acquisition Site" : "collection_site", "Age At Acquisition (Years)" : "age_in_years_at_collection", "TNM Stage" : "stage", "Histological Grade" : "grade" }
    HCMI_sample_columns = hcmiDf.loc[:, sample_columns].rename(columns=sample_name_map)
    HCMI_sample_columns["patient_id"] = HCMI_sample_columns['model_id'].apply(generate_patient_id)
    HCMI_sample_columns["sample_id"] = HCMI_sample_columns['model_id'].apply(generate_sample_id)
    HCMI_sample_columns["diagnosis"] = HCMI_sample_columns["Clinical Tumor Diagnosis"] + " " + HCMI_sample_columns["Histological Subtype"]
    HCMI_sample_columns = HCMI_sample_columns.drop(["Clinical Tumor Diagnosis","Histological Subtype"], axis=1)
    HCMI_sample_columns["staging_system"] = "TNM staging system"
    HCMI_sample_columns["grading_system"] = "AJCC"
    HCMI_sample_columns['sharable'] = 'yes'
    HCMI_sample_columns['treated_prior_to_collection'] = clinicalDf['prior_treatment']
    #import pdb; pdb.set_trace()
    indexOfPriorTreated = HCMI_sample_columns['treated_prior_to_collection'].str.lower() == "yes"
    notIndexOfPriorTreated = HCMI_sample_columns['treated_prior_to_collection'].str.lower() != "yes"
    HCMI_sample_columns.loc[indexOfPriorTreated, "treatment_naive_at_collection"] = "no"
    HCMI_sample_columns.loc[notIndexOfPriorTreated, "treatment_naive_at_collection"] = "not collected" 
    accepted_tumour_types = ["Metastasis", "Primary", "Recurrent"]
    HCMI_sample_columns = HCMI_sample_columns[HCMI_sample_columns.loc[:,'tumour_type'].isin(accepted_tumour_types)]
    align_and_append_sheet(sample_template, HCMI_sample_columns).to_csv("HCMI_metadata-patient_sample.tsv", sep='\t', index=False)

def generate_other_models_sheet(hcmiDf, other_model_template):
    other_model_columns = ["Name", "Link To Distributor", "Type"]
    other_model_columns_rename_map = {"Name":"name", "Link To Distributor": "supplier", "Type":"type"}
    HCMI_other_models_sheet = hcmiDf.loc[:, other_model_columns].rename(columns=other_model_columns_rename_map)
    HCMI_other_models_sheet['supplier'] = "ATCC:" + HCMI_other_models_sheet['supplier']
    HCMI_other_models_sheet['model_id'] = HCMI_other_models_sheet['name']
    align_and_append_sheet(other_model_template, HCMI_other_models_sheet)\
        .dropna(axis=1, how='all').to_csv("HCMI_metadata-other_models.tsv", sep='\t', index=False)

def generate_sharing(hcmiDf, sharing_template):
    sharing_columns = ["Name", "Link to Model Details", "Licensing Required For Commercial Use"]
    sharing_columns_rename_map = {"Name":"model_id", "Licensing Required For Commercial Use": "accessibility" }
    contact_email = "ocg@mail.nih.gov"
    HCMI_sharing_sheet = hcmiDf.loc[:, sharing_columns].rename(columns=sharing_columns_rename_map)
    HCMI_sharing_sheet['email'] = contact_email
    HCMI_sharing_sheet['accessibility'] = "industry and academia"
    HCMI_sharing_sheet['database_url'] = f'https://hcmi-searchable-catalog.nci.nih.gov/model/{model_id}'
    align_and_append_sheet(sharing_template, HCMI_sharing_sheet).to_csv("HCMI_metadata-sharing.tsv", sep='\t', index=False)

hcmiDf = pd.read_csv("model-table.tsv", sep='\t')
clinicalDf = pd.read_csv("clinical.tsv", sep='\t')
patient_template = pd.read_csv("metadata_template-patient.tsv", sep='\t').dropna(axis=0, how='all')
sample_template = pd.read_csv("metadata_template-patient_sample.tsv", sep='\t').dropna(axis=0, how='all')
other_models_template = pd.read_csv("metadata_template-other_models.tsv", sep='\t').dropna(axis=0, how='all')
model_validation = pd.read_csv("metadata_template-model_validation.tsv", sep='\t').dropna(axis=0, how='all')
sharing_template = pd.read_csv("metadata_template-sharing.tsv", sep='\t').dropna(axis=0, how='all')
generate_patient_sheet(hcmiDf, patient_template) 
generate_sample_sheet(hcmiDf, sample_template, clinicalDf)
generate_other_models_sheet(hcmiDf, other_models_template)
generate_sharing(hcmiDf, sharing_template)
