import pandas as pd
from numpy import where
from re import sub
from os import listdir, getcwd, rename, makedirs
from os.path import isfile, join, isdir, exists


def read_tsv(path):
    return pd.read_csv(path, delimiter='\t', na_values='', skiprows=(1,2,3,4))

def write_tsv(path, df):
    df.to_csv(path , sep='\t', index=False)

def read_data(path):
    patient_data = pd.read_excel(io = path, sheet_name=0, header=0, skiprows=[1],keep_default_na=False)
    clinical_data = pd.read_excel(io = path, sheet_name=1, header=0, skiprows=[1],keep_default_na=False)
    treatment_data = pd.read_excel(io = path, sheet_name=2, header=0, skiprows=[1],keep_default_na=False)
    model_data = pd.read_excel(io = path, sheet_name=3, header=0, skiprows=[1],keep_default_na=False)
    return patient_data, clinical_data, treatment_data, model_data

def transform_patient_data(patient):
    patient['ethnicity_assessment_method'] = ""
    out_patient = patient[['Patient ID *', 'Gender *', 'History of Smoking', 'Ethnicity *', 'ethnicity_assessment_method']]
    out_patient['initial_diagnosis'] = ""
    out_patient['age_at_initial_diagnosis'] = ""
    SH_rep = {'Never': 'Smoking History:non-smoker', 'Former': 'Smoking History:ex-smoker', 'Current':'Smoking History:current-smoker'}
    out_patient.replace(SH_rep, inplace=True)
    mapper = {'Patient ID *': 'patient_id', 'Gender *': 'sex', 'History of Smoking': 'history', 'Ethnicity *':'ethnicity'}
    return out_patient.rename(mapper, axis=1)

def merge_stage(value):
    string = ''
    if value[0] != '':
        string+=value[0]
    if value[1] !='':
        if string!='':
            string+=','
        string+=value[1]
    if value[2] !='':
        if string!='':
            string+=','
        string+=value[2]
    return string

def extract_grade(row):
    if row=='':
        return 'Not provided'
    else:
        if row =='High Grade' or row =='Low Grade':
            return row
        else:
            row = sub(r" (.*)", "", row)
            return sub(r"G","Grade ", row)

def transform_patient_sample_data(path, clinical, model):
    ## Patient sample
    # patient_id	sample_id	collection_date	collection_event	months_since_collection_1	age_in_years_at_collection	diagnosis	tumour_type	primary_site	
    # collection_site	stage	staging_system	grade	grading_system	virology_status	sharable	treatment_naive_at_collection	
    # treated_at_collection	treated_prior_to_collection	model_id

    clinical['age'] = clinical['Event Age *']
    model['age'] = model['Age at Collection *'] 
    merged = pd.merge(clinical, model, on=['Patient ID *', 'age'])
    dups = merged[merged.duplicated(subset=['Patient ID *', 'age'], keep=False)]
    print("Number of duplicate entries: %d" %len(dups))
    #write_tsv(join(path, 'Dups.tsv'), merged)
    merged['tumor_type'] = merged['Event Setting *'] ## what about normal/benign
    merged['primary_site'] = merged['Disease Group *']
    merged['collection_site'] = merged['Specimen Site'] ## should replace Pleural effusion to lung?
    merged['stage'] = merged[['pT','pN','pM' ]].apply(merge_stage, axis=1)
    merged['staging_system'] = where(merged['stage'].eq(''), '', 'TNM staging system')
    merged['grade'] = merged['Tumor Grade'].apply(extract_grade)
    merged['grading_system'] = 'Not provided'
    out_patient_sample = merged
    return out_patient_sample



input_path = "E:/EBI/Work/PDXNet/PDXPortal_Breast_Template_to_HCI-054.xlsx"
meta_path = "E:/EBI/Work/PDCM_data/active_templates/metadata"

patient, clinical, treatment, model = read_data(input_path)

m_model, m_patient, m_p_sample, m_pdx, m_share = read_tsv(join(meta_path, 'metadata_template-model_validation.tsv')), read_tsv(join(meta_path, 'metadata_template-patient.tsv')), \
                                                 read_tsv(join(meta_path, 'metadata_template-patient_sample.tsv')), read_tsv(join(meta_path, 'metadata_template-pdx_models.tsv')), read_tsv(join(meta_path, 'metadata_template-sharing.tsv'))
                                                 
m_model, m_patient, m_p_sample, m_pdx, m_share  = m_model.drop(['Field'], axis=1), m_patient.drop(['Field'], axis=1), m_p_sample.drop(['Field'], axis=1), m_pdx.drop(['Field'], axis=1), m_share.drop(['Field'], axis=1)


    
#def transform_pdx_data(pdx):
    ## PDX data
    # model_id	host_strain_name	host_strain_nomenclature	engraftment_site	engraftment_type	sample_type	sample_state	passage_number	publications
    

#merged = pd.merge(patient[['Patient ID *', 'Gender *', 'Ethnicity *', 'History of Smoking']], clinical[['Patient ID *', 'Event Age *']], on='Patient ID *', how='left').sort_values('Event Age *', ascending=True).drop_duplicates('Patient ID *').reset_index(drop=True)
#merged['ethnicity_assessment_method'] = ""