#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import sys
import math
import os

def align_and_append_sheet(template, other_sheet):
    alignment = template.align(other_sheet, join="outer")[1]
    return template.append(alignment).dropna(axis=0, how='all').dropna(axis=1, how='all')

def numToMon(num):
    month_str = ""
    if not math.isnan(num):
        num_int = int(num)
        months = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec" ]
        month_str = months[num_int - 1]
    return month_str

def transformCellPassport_patient(cellPassportDf, patient_template):
    patient_columns = ["patient_id", "gender", "smoking_status", "ethnicity"]
    patient_column_name_changes = {"gender":"sex", "smoking_status":"history"}
    cell_passport_patient_sheet = cellPassportDf.loc[:, patient_columns].rename(columns = patient_column_name_changes)
    save_to_file = "pdcm_format/CMP_metadata-patient_sample.tsv"
    print(f'Saving {save_to_file}')
    align_and_append_sheet(patient_template, cell_passport_patient_sheet).to_csv(save_to_file, sep='\t', index=False)


def transformCellPassport_sample(cellPassportDf, sample_template):
    sample_columns = ["patient_id", "sample_id", "model_id", "sampling_year", "sampling_month", "age_at_sampling", "cancer_type_detail", "tissue_status", "tissue", "sample_site", "tnm_t", "tnm_n", "tnm_m", "tumour_grade", "sample_treatment"]
    cell_passport_columns = cellPassportDf.loc[:, sample_columns]
    columns_To_Drop = ["sampling_year", "sampling_month", "tnm_t", "tnm_n", "tnm_m", "sample_treatment"]
    sample_rename_mapping = { "age_at_sampling" : "age_in_years_at_collection", "cancer_type_detail":"diagnosis", "tissue_status":"tumour_type", "tissue" : "primary_site", "sample_site":"collection_site", "tumour_grade":  "grade" }
    cell_passport_sample_sheet = cell_passport_columns.drop(columns=columns_To_Drop).rename(columns=sample_rename_mapping)
    #Format collection date
    cell_passport_columns["sampling_month"] = cell_passport_columns["sampling_month"].apply(numToMon)
    cell_passport_sample_sheet['collection_date'] = cell_passport_columns["sampling_month"] + " " + cell_passport_columns["sampling_year"].astype(str)
    cell_passport_sample_sheet['collection_date'] = cell_passport_sample_sheet['collection_date'].str.replace(" nan", "")
    #Cat TNM staging system
    cell_passport_sample_sheet['stage'] =  cell_passport_columns["tnm_t"] + "," + cell_passport_columns["tnm_n"] + "," +  cell_passport_columns["tnm_m"]
    #Hard coded values
    cell_passport_sample_sheet["staging_system"] = "TNM staging system"
    cell_passport_sample_sheet["sharable"] = "yes"
    #Need to proccess sample_treatment as well still
    save_to_file = "pdcm_format/CMP_metadata-sample.tsv"
    print(f'Saving {save_to_file}')
    align_and_append_sheet(sample_template,cell_passport_sample_sheet).to_csv(save_to_file, sep='\t', index=False)


def transformedCellPassport_cell_models(cellPassportDf, cell_models_template):
    other_models_columns = ["model_id", "sample_id", "parent_id", "model_name", "model_type", "growth_properties",                            "model_comments", "model_relations_comment", "suppliers", "COSMIC_ID", "BROAD_ID", "RRID",                             "CCLE_ID", "pmed"]
    cp_other_models_columns =  cellPassportDf.loc[:, other_models_columns]
    other_models_rename_map = {"model_name":"name", "model_type":"type", "sample_id": "origin_patient_sample_id",                             "suppliers":"supplier", "pmed":"publications" }
    cp_cell_models_sheet = cp_other_models_columns.rename(columns=other_models_rename_map)
    cp_cell_models_sheet['comments'] = cp_other_models_columns[["model_comments", "model_relations_comment"]].fillna('').agg(' '.join, axis=1).str.replace("\n"," ")
    cp_other_models_columns['prepended_cosmic_id'] = cp_other_models_columns['COSMIC_ID'].astype(str).apply(lambda x: if_str_not_na_prepend_str(x, "COSMIC"))
    cp_other_models_columns['prepended_CCLE_id'] = cp_other_models_columns['CCLE_ID'].astype(str).apply(lambda x: if_str_not_na_prepend_str(x, "CCLE_Name:"))
    cp_cell_models_sheet['external_ids'] = cp_other_models_columns[['prepended_cosmic_id', 'BROAD_ID', 'RRID', "prepended_CCLE_id"]]        .fillna('')        .astype(str)        .agg(','.join, axis=1)
    cp_cell_models_sheet['external_ids'] = cp_cell_models_sheet['external_ids'].str.replace(r'^,+', '', regex=True)
    cp_cell_models_sheet['external_ids'] = cp_cell_models_sheet['external_ids'].str.replace(r',+$', '', regex=True)
    cp_cell_models_sheet.loc[cp_cell_models_sheet["parent_id"].isna() ,"origin_patient_sample_id"] = ""
    cp_cell_models_sheet.loc[cp_cell_models_sheet["parent_id"].isna() ,"origin_patient_sample_id"] = ""
    columns_to_drop = ["model_comments", "model_relations_comment", "COSMIC_ID", "BROAD_ID", "RRID", "CCLE_ID"]
    cp_cleaned_cell_model_sheet = cp_cell_models_sheet.drop(columns= columns_to_drop)
    save_to_file = "pdcm_format/CMP_metadata-cell_model.tsv"
    print(f'Saving {save_to_file}')
    align_and_append_sheet(cell_models_template, cp_cleaned_cell_model_sheet).to_csv(save_to_file, sep='\t', index=False)

    pass

def transformCellPassport_sharing(cellPassportDf, sharing_template):
    CMP_sharing_sheet = pd.DataFrame()
    CMP_sharing_sheet['model_id'] = cellPassportDf['model_id']
    CMP_sharing_sheet['accessibility'] = 'academia and industry'
    sharing_filename = "pdcm_format/CMP_metadata-sharing.tsv"
    print(f'Saving {sharing_filename}')
    align_and_append_sheet(sharing_template, CMP_sharing_sheet).to_csv(sharing_filename, sep='\t', index=False)


def if_str_not_na_prepend_str(str_to_check, str_to_prepend):
    built_str = ''
    str_to_check_arr = str_to_check.split(";")
    space = " "
    count = 0
    for str_to_check in str_to_check_arr:
        if str_to_check != 'nan' and str_to_check_arr != '' and str_to_check:
            built_str += f'{str_to_prepend}{str_to_check}{space*count}'
    return built_str


def getTemplates():
    patient_template = pd.read_csv("metadata_template-patient.tsv", sep='\t').dropna(axis=0, how='all')
    sample_template = pd.read_csv("metadata_template-sample.tsv", sep='\t').dropna(axis=0, how='all')
    cell_models_template = pd.read_csv("metadata_template-cell_model.tsv", sep='\t').dropna(axis=0,
                                                                                                                 how='all')
    model_validation_template = pd.read_csv("metadata_template-model_validation.tsv", sep='\t').dropna(axis=0,
                                                                                                       how='all')
    sharing_template = pd.read_csv("metadata_template-sharing.tsv", sep='\t').dropna(axis=0, how='all')
    return patient_template, sample_template, cell_models_template, model_validation_template, sharing_template


def filterMainDataFrame(cellPassportDf):
    tumour_type_filter = ["Metastasis", "Unknown", "Tumour"]
    filterIndex = cellPassportDf['HCMI'].isna() & cellPassportDf['tissue_status'].isin(tumour_type_filter) & \
                  cellPassportDf['suppliers'].notna() & cellPassportDf['suppliers'].str.match("(?!Unknown:Unknown)")
    return cellPassportDf[filterIndex].copy()

def transformCellPassport_model_validation(filteredCpDf, model_validation_template):
    return False

def createCellPassportPdxFinderTemplates():
    cellPassportDf = pd.read_csv("model_list_20220810.csv", engine="c", sep=',').dropna(axis=0, how='all')
    if not os.path.exists('pdcm_format'):
        os.makedirs('pdcm_format')
    filteredCpDf = filterMainDataFrame(cellPassportDf)
    patient_template, sample_template, cell_models_template, model_validation_template, sharing_template = getTemplates()
    transformCellPassport_patient(filteredCpDf, patient_template)
    transformCellPassport_sample(filteredCpDf, sample_template)
    transformedCellPassport_cell_models(filteredCpDf, cell_models_template)
    transformCellPassport_model_validation(filteredCpDf, model_validation_template)
    transformCellPassport_sharing(filteredCpDf, sharing_template)

createCellPassportPdxFinderTemplates()


