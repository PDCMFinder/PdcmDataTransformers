import pandas as pd


def compare_pdx(old_path, new_path):
    old = pd.read_csv(old_path, sep='\t', na_values="", low_memory=False, skiprows=[1,2,3,4]).drop('Field',axis=1)
    old.model_id = old.model_id.replace(to_replace='PHLC0', value='PHLC',regex=True)
    old_model_id = old.model_id.unique()

    new = pd.read_csv(new_path, sep='\t', na_values="", low_memory=False)
    new.model_id = new.model_id.replace(to_replace=' ', value='',regex=True)
    new_model_id = new.model_id.unique()

    subset = set(old_model_id) - set(new_model_id)
    subset = old[old.model_id.isin(list(subset))]

    subset = pd.concat([subset,new])
    return subset
old_path = "~/pdx/pdxfinder-data/data/UPDOG/PMLB/PMLB_metadata-pdx_model.tsv"
new_path = "~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/PMLB_metadata-pdx_model.tsv"
subset = compare_pdx(old_path, new_path)
subset.to_csv('~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/metadata_updated/PMLB_metadata-pdx_model.tsv', sep='\t', index=False)

old_path = "~/pdx/pdxfinder-data/data/UPDOG/PMLB/PMLB_metadata-model_validation.tsv"
new_path = "~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/PMLB_metadata-model_validation.tsv"
subset = compare_pdx(old_path, new_path)
subset.to_csv('~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/metadata_updated/PMLB_metadata-model_validation.tsv', sep='\t', index=False)

old_path = "~/pdx/pdxfinder-data/data/UPDOG/PMLB/PMLB_metadata-patient_sample.tsv"
new_path = "~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/PMLB_metadata-patient_sample.tsv"
subset = compare_pdx(old_path, new_path)
subset.to_csv('~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/metadata_updated/PMLB_metadata-patient_sample.tsv', sep='\t', index=False)

old_path = "~/pdx/pdxfinder-data/data/UPDOG/PMLB/PMLB_metadata-sharing.tsv"
new_path = "~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/PMLB_metadata-sharing.tsv"
subset = compare_pdx(old_path, new_path)
subset.to_csv('~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/metadata_updated/PMLB_metadata-sharing.tsv', sep='\t', index=False)

def compare_patient_sample(old_path, new_path):
    old = pd.read_csv(old_path, sep='\t', na_values="", low_memory=False, skiprows=[1,2,3,4]).drop('Field',axis=1)
    old_patient_id = old.patient_id.unique()

    new = pd.read_csv(new_path, sep='\t', na_values="", low_memory=False)
    new_patient_id = new.patient_id.unique()

    subset = set(old_patient_id) - set(new_patient_id)
    subset = old[old.patient_id.isin(list(subset))]

    subset = pd.concat([subset,new])
    return subset

old_path = "~/pdx/pdxfinder-data/data/UPDOG/PMLB/PMLB_metadata-patient.tsv"
new_path = "~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/PMLB_metadata-patient.tsv"
compare_patient_sample(old_path, new_path)
subset.to_csv('~/pdx/update-data/data-repo/PMLB_April_2022/PMLB/metadata_updated/PMLB_metadata-patient.tsv', sep='\t', index=False)



