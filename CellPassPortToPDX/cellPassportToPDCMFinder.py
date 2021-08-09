#!/usr/bin/env python
# coding: utf-8

# In[134]:


import pandas as pd
import sys
import math


# In[135]:


from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[172]:


cellPassportDf = pd.read_excel("model_list_20210310.xlsx", engine="openpyxl")


# In[173]:


patient_template = pd.read_csv("metadata_template-patient.tsv", sep='\t').dropna(axis=0, how='all')
sample_template = pd.read_csv("metadata_template-sample.tsv", sep='\t').dropna(axis=0, how='all')


# In[174]:


def align_and_append_sheet(template, other_sheet):
    alignment = template.align(other_sheet, join="outer")[1]
    return template.append(alignment).dropna(axis=0, how='all')


# In[175]:


def transformCellPassport_patient(cellPassportDf):
    patient_columns = ["patient_id", "gender", "smoking_status", "ethnicity"]
    # No direct value changes are needed for the patient sheet
    patient_column_name_changes = {"gender":"sex", "smoking_status":"history"}
    cell_passport_patient_sheet = cellPassportDf.loc[:, patient_columns].rename(columns = patient_column_name_changes)
    align_and_append_sheet(patient_template, cell_passport_patient_sheet).to_csv("HCMI_metadata-patient.tsv", sep='\t', index=False)


# In[176]:


def numToMon(num):
    month_str = ""
    if not math.isnan(num):
        num_int = int(num)
        months = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec" ]
        month_str = months[num_int - 1]
    return month_str
        


# In[177]:


sample_columns = ["patient_id", "sample_id", "model_id", "sampling_year", "sampling_month", "age_at_sampling", "cancer_type_detail", "tissue_status", "tissue", "sample_site", "tnm_t", "tnm_n", "tnm_m", "tumour_grade", "sample_treatment"]
cell_passport_columns = cellPassportDf.loc[:, sample_columns] 
cell_passport_sample_sheet = pd.DataFrame()

columns_To_Drop = ["sampling_year", "sampling_month", "tnm_t", "tnm_n", "tnm_m", "sample_treatment"]
sample_rename_mapping = { "age_at_sampling" : "age_in_years_at_collection", "cancer_type_detail":"diagnosis", "tissue_status":"tumour_type", "tissue" : "primary_site", "sample_site":"collection_site", "tumour_grade":  "grade" }
cell_passport_sample_sheet = cell_passport_columns.drop(columns=columns_To_Drop).rename(columns=sample_rename_mapping)
#Format collection date
cell_passport_columns["sampling_month"] = cell_passport_columns["sampling_month"].apply(numToMon)


# In[178]:


cell_passport_sample_sheet['collection_date'] = cell_passport_columns["sampling_month"] + " " + cell_passport_columns["sampling_year"].astype(str)
cell_passport_sample_sheet['collection_date'] = cell_passport_sample_sheet['collection_date'].str.replace(" nan", "")
#Cat TNM staging system
cell_passport_sample_sheet['stage'] =  cell_passport_columns["tnm_t"] + "," + cell_passport_columns["tnm_n"] + "," +  cell_passport_columns["tnm_m"]
#Hard coded values
cell_passport_sample_sheet["staging_system"] = "TNM staging system"
cell_passport_sample_sheet["sharable"] = "yes"


# In[179]:


tumour_type_filter = ["Metastasis", "Unknown", "Tumour"]
filtered_passport = cell_passport_sample_sheet[cell_passport_sample_sheet["tumour_type"].isin(tumour_type_filter)]

#Need to proccess sample_treatment as well still
align_and_append_sheet(sample_template,filtered_passport).to_csv("HCMI_metadata-sample.tsv", sep='\t', index=False)


# In[181]:


used_columns = ["gender", "smoking_status", "ethnicity",  "sample_id", "model_id", "sampling_year", "sampling_month", "age_at_sampling", "cancer_type_detail", "tissue_status", "tissue", "sample_site", "tnm_t", "tnm_n", "tnm_m", "tumour_grade", "sample_treatment"]
cellPassportDf.drop(columns=used_columns).to_csv("leftovers.tsv", sep='\t', index=False)


# In[ ]:




