import pandas as pd
from numpy import where
from re import sub
#from os import listdir, getcwd, rename, makedirs
from os.path import join


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

def generate_sample_id(row):
    if row[0] == "":
        row[0]="U" ## for missing event type
    return 'SID-'+str(row[0][0])+'-'+str(row[1])

def extract_engraftment_type(row):
    value=''
    if row[1]!= '' or row[0]!='':
        if 'mammary fat' in row[1]:
            if 'Breast' in row[0]:
                value='orthotopic'    
            else: 
                value='heterotopic'
    else:
        value='Not provided'
    return value

def extract_parent_pdx_model(parent_model):
    if len(parent_model == 1):
        m_id = str(parent_model.engraftment_type.str.replace(' PDX', ''))
    else:
        m_id = ''
        print('Check manually: Pateint_id: %s' %parent_model.pateint_id[0])    
    return m_id

def generate_passage_number(row, pdx):
    ## Still needs to be fixed. Buggy.
    if row[3] != 0:
        p_id = row[0]
        subset = pdx[pdx.patient_id == p_id]
        if len(subset) > 1:
            if any(subset.engraftment_type.str.contains('PDX')):
                parent_model = subset[subset.engraftment_type.str.contains('PDX')]
                m_id = extract_parent_pdx_model(parent_model)
                if m_id != '':
                    if m_id in subset.model_id:
                        return 0
                        
                    
        else:
            if any(subset.engraftment_type.str.contains('PDX')):
                parent_model = subset[subset.engraftment_type.str.contains('PDX')]
                m_id = extract_parent_pdx_model(parent_model)            
    else:
        return row[3]


def transform_patient_data(patient):
    patient['ethnicity_assessment_method'] = ""
    out_patient = patient[['Patient ID *', 'Gender *', 'History of Smoking', 'Ethnicity *', 'ethnicity_assessment_method']]
    out_patient['initial_diagnosis'] = ""
    out_patient['age_at_initial_diagnosis'] = ""
    SH_rep = {'Never': 'Smoking History:non-smoker', 'Former': 'Smoking History:ex-smoker', 'Current':'Smoking History:current-smoker'}
    out_patient.replace(SH_rep, inplace=True)
    mapper = {'Patient ID *': 'patient_id', 'Gender *': 'sex', 'History of Smoking': 'history', 'Ethnicity *':'ethnicity'}
    return out_patient.rename(mapper, axis=1)
    
def transform_patient_sample_data(path, clinical, model):
    ## Patient sample
    # patient_id	sample_id	collection_date	collection_event	months_since_collection_1	age_in_years_at_collection	diagnosis	tumour_type	primary_site	
    # collection_site	stage	staging_system	grade	grading_system	virology_status	sharable	treatment_naive_at_collection	
    # treated_at_collection	treated_prior_to_collection	model_id

    clinical['age'] = clinical['Event Age *']
    model['age'] = model['Age at Collection *'] 
    merged = pd.merge(clinical[['Patient ID *', 'Event Type *', 'age', 'Event Setting *', 'pT', 'pN', 'pM', 'Tumor Grade', 'Treatment Naïve']], model[['Patient ID *', 'age', 'Model ID *', 'Disease Group *', 'Specimen Site']], on=['Patient ID *', 'age'])
    dups = merged[merged.duplicated(subset=['Patient ID *', 'age'], keep=False)]
    print("Number of duplicate entries: %d" %len(dups))
    write_tsv(join(path, 'Dups.tsv'), merged)
    merged['tumour_type'] = merged['Event Setting *'].str.replace(r'Distant | Second ', '', regex=True)\
                                                     .str.replace(r'Metastasis', 'metastastic', regex=True)\
                                                     .str.replace(r'Recurrence', 'recurrent', regex=True)\
                                                         .str.replace(r'Primary', 'primary', regex=True)## what about normal/benign
    
    merged['primary_site'] = merged['Disease Group *']
    merged['collection_site'] = merged['Specimen Site'] ## should replace Pleural effusion to lung?
    missing_cs = merged['collection_site'].str.contains("''|PDX")
    print("Number of entries with missing collection site: %d" %len(missing_cs))
    write_tsv(join(path, 'missing_cs.tsv'), merged)
    merged['stage'] = merged[['pT','pN','pM' ]].apply(merge_stage, axis=1)
    merged['staging_system'] = where(merged['stage'].eq(''), '', 'TNM staging system')
    merged['grade'] = merged['Tumor Grade'].apply(extract_grade)
    merged['grading_system'] = 'Not provided' ## Ask
    merged['virology_status'] = 'Not provided' ## Ask
    merged['sharable'] = 'Not provided' ## Ask
    merged['treatment_naive_at_collection'] = merged['Treatment Naïve'].replace('', 'Not provided')
    merged['treated_at_collection'] = "Not provided" ## Ask
    merged['treated_prior_to_collection'] = 'Not provided' ## Ask
    merged['sample_id'] = merged[['Event Type *', 'Model ID *']].apply(generate_sample_id, axis=1)
    mapper ={'Patient ID *': 'patient_id', 'Model ID *': 'model_id', 'age': 'age_in_years_at_collection'}
    out_patient_sample = merged.rename(mapper, axis = 1)
    out_patient_sample['diagnosis'] = out_patient_sample['Disease Group *'] + ' cancer'
    out_patient_sample['collection_date'] = ""
    out_patient_sample['collection_event'] = ""
    out_patient_sample['months_since_collection_1'] = ""
    out_patient_sample = out_patient_sample[["patient_id", "sample_id", "collection_date", "collection_event", "months_since_collection_1", "age_in_years_at_collection", "diagnosis", "tumour_type", "primary_site", "collection_site", "stage", "staging_system", "grade", "grading_system", "virology_status", "sharable", "treatment_naive_at_collection", "treated_at_collection", "treated_prior_to_collection", "model_id"]]
    return out_patient_sample

def transform_pdx_data(pdx):
    ## PDX data
    # model_id	host_strain_name	host_strain_nomenclature	engraftment_site	engraftment_type	sample_type	sample_state	
    # passage_number	publications
    mapper ={'Patient ID *': 'patient_id', 'Model ID *': 'model_id', 'Host Strain':'host_strain_name', 
             'Engraftment Site':'engraftment_site', 'Specimen Condition': 'sample_state', 
                 'Specimen Site': 'engraftment_type', 'Engraftment Material': 'sample_type'}
    pdx = pdx.rename(mapper, axis=1)
    pdx['sample_state'] = pdx['sample_state'].replace("", 'Not provided')
    pdx['passage_number'] = ''
    pdx.passage_number[pdx.engraftment_type.str.contains(' PDX') != True] = 0
    subset = pdx[['patient_id', 'model_id', 'engraftment_type', 'passage_number']]
    #pdx['passage_number'] = [generate_passage_number(row[1], subset) for row in subset.iterrows()]
    pdx['engraftment_type'] = pdx[['engraftment_type', 'engraftment_site']].apply(extract_engraftment_type, axis=1)
    pdx['host_strain_nomenclature'] = ""
    pdx['publications'] = ""
    pdx = pdx[['model_id', 'host_strain_name', 'host_strain_nomenclature', 'engraftment_site', 'engraftment_type', 'sample_type', 'sample_state', 'passage_number', 'publications']]
    return pdx

def transform_sharable(pdx):
    ## model_id	accessibility	europdx_access_modality	email	name	form_url	database_url
    mapper ={'Patient ID *': 'patient_id', 'Model ID *': 'model_id', 'Public Status *':'accessibility'}
    pdx = pdx.rename(mapper, axis=1)
    sharing = pdx[['model_id', 'accessibility']]
    sharing[['europdx_access_modality', 'email', 'name', 'form_url', 'database_url']] = ''
    return sharing

def transform_model_validation(model):
    mv = model['Model ID *']
    mv[['validation_technique','description','passages_tested','validation_host_strain_nomenclature']] = ""
    return mv
 
def flatten_cytogenetics(temp, symbol):
    temp['symbol'] = symbol
    temp[['essential_or_additional_marker','platform_id']] = ''
    return temp

def transform_cytogenetics(clinical, patient_sample):
    ## sample_id	symbol	marker_status	essential_or_additional_marker	platform_id
    cyto = clinical[['Patient ID *', 'ER Status', 'PR Status', 'HER2 Status', 'KI67 Status', 'PAM50 Subtype', 'Event Age *']]
    mapper = {'Patient ID *':'patient_id', 'Event Age *': 'age', 'ER Status':'er_status', 'PR Status':'pr_status', 'HER2 Status':'her2_status', 'KI67 Status':'ki67_status', 'PAM50 Subtype':'pam50_subtype'}
    cyto = cyto.rename(mapper, axis=1)
    patient_sample['age'] = patient_sample['age_in_years_at_collection']
    cyto = pd.merge(cyto, patient_sample[['patient_id','model_id', 'sample_id', 'age']], on=['patient_id', 'age'])
    cyto = cyto.drop([16,17],axis=0).reset_index(drop=True)
    cyto = cyto.drop(['age', 'model_id'],axis=1)
    #out_cyto['sample_id'] = cyto['sample_id']
    ## ER
    out_cyto = flatten_cytogenetics(cyto[['sample_id', 'er_status']], 'ESR1').rename({'er_status':'marker_status'}, axis=1)
    ##HER 
    out_cyto = pd.concat([out_cyto, flatten_cytogenetics(cyto[['sample_id', 'her2_status']], 'ERBB2').rename({'her2_status':'marker_status'}, axis=1)]).reset_index(drop=True)
    out_cyto = pd.concat([out_cyto, flatten_cytogenetics(cyto[['sample_id', 'pr_status']], 'PGR').rename({'pr_status':'marker_status'}, axis=1)]).reset_index(drop=True)
    out_cyto.marker_status = out_cyto.marker_status.replace(r'^\s*$', 'Not provided', regex=True)
    out_cyto['platform_id'] = 'cytogenetics_immunohistochemistry'
    #out_cyto = pd.concat([out_cyto, flatten_cytogenetics(cyto[['sample_id', 'ki67_status']], 'ERBB2').rename({'her2_status':'marker_status'}, axis=1)]).reset_index(drop=True)
    return out_cyto

    
    
def generate_tsv(input_path):
    patient, clinical, treatment, model = read_data(input_path)
    model = model.iloc[0:18,]
    write_tsv(join(output_path, 'PDXNet_metadata_patient.tsv'), transform_patient_data(patient))
    ps_models = transform_patient_sample_data(output_path, clinical, model)
    write_tsv(join(output_path, 'PDXNet_metadata_patient_sample.tsv'), ps_models)
    write_tsv(join(output_path, 'PDXNet_metadata_pdx_model.tsv'), transform_pdx_data(model))
    #write_tsv(join(output_path, 'PDXNet_metadata_model_validation.tsv'), transform_model_validation(model))
    write_tsv(join(output_path, 'PDXNet_metadata_sharing.tsv'), transform_sharable(model))
    
    
    ## Before cyto:
    patient_sample = ps_models.drop([16,17],axis=0).reset_index(drop=True)
    write_tsv(join(output_path, 'PDXNet_cytogenetics.tsv'), transform_cytogenetics(clinical, patient_sample))
    
    
input_path = "E:/EBI/Work/PDXNet/PDXPortal_Breast_Template_to_HCI-054.xlsx"
output_path = "E:/EBI/Work/PDXNet/tsv"
meta_path = "E:/EBI/Work/PDCM_data/active_templates/metadata"
generate_tsv(input_path)

#m_model, m_patient, m_p_sample, m_pdx, m_share = read_tsv(join(meta_path, 'metadata_template-model_validation.tsv')), read_tsv(join(meta_path, 'metadata_template-patient.tsv')), \
#                                                 read_tsv(join(meta_path, 'metadata_template-patient_sample.tsv')), read_tsv(join(meta_path, 'metadata_template-pdx_models.tsv')), read_tsv(join(meta_path, 'metadata_template-sharing.tsv'))
                                                 
#m_model, m_patient, m_p_sample, m_pdx, m_share  = m_model.drop(['Field'], axis=1), m_patient.drop(['Field'], axis=1), m_p_sample.drop(['Field'], axis=1), m_pdx.drop(['Field'], axis=1), m_share.drop(['Field'], axis=1)
#merged = pd.merge(patient[['Patient ID *', 'Gender *', 'Ethnicity *', 'History of Smoking']], clinical[['Patient ID *', 'Event Age *']], on='Patient ID *', how='left').sort_values('Event Age *', ascending=True).drop_duplicates('Patient ID *').reset_index(drop=True)
#merged['ethnicity_assessment_method'] = ""