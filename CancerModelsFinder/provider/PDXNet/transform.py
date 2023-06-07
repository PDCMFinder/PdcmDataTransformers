from CancerModelsFinder.provider.PDXNet.resources import *
from CancerModelsFinder.classes import transform

from os.path import join, exists
import pandas as pd


class PDXNet_data_transformations():
    file_paths = file_paths

    def __init__(self):
        self.flag = None
        self.PDXNet_models = None
        self.ref_data = pd.read_csv(self.file_paths['JAX_endpoint_reference_file'], sep='\t')
        self.transformations = transform.data_transformations()

    def process_models(self):
        if exists(self.file_paths['model']) and exists(self.file_paths['patient']) and exists(self.file_paths['sample']) \
                and exists(self.file_paths['PDMR_mol']) and exists(self.file_paths['PDTC_mol']):
            model = pd.read_csv(self.file_paths['model']).drop('Unnamed: 0', axis=1)
            patient = pd.read_csv(self.file_paths['patient']).drop('Unnamed: 0', axis=1)
            sample = pd.read_csv(self.file_paths['sample']).drop('Unnamed: 0', axis=1)
            PDMR_mol = pd.read_csv(self.file_paths['PDMR_mol']).drop('Unnamed: 0', axis=1)
            PDTC_mol = pd.read_csv(self.file_paths['PDTC_mol']).drop('Unnamed: 0', axis=1)
            self.PDXNet_models = [model, patient, sample, PDMR_mol, PDTC_mol]

    def models_in_CM(self):
        models_in_cancerModels = pd.DataFrame()
        for provider in ['HCI-BCM', 'MDAnderson', 'WISTAR-MDAnderson-Penn', 'WUSTL']:
            models = pd.read_csv(join(self.file_paths['source'], provider, provider + '_metadata-pdx_model.tsv'),
                                 sep='\t')
            models = self.transformations.drop_template_rows_columns(models)
            models['Provider'] = provider
            models_in_cancerModels = pd.concat([models_in_cancerModels, models])
        self.model_list_in_CM = models_in_cancerModels

    def compare_models(self):
        if self.PDXNet_models is not None:
            self.models_in_CM()
            old_models = self.model_list_in_CM['model_id'].to_list()
            PDXNet_model_list = self.PDXNet_models[1]['ContributorPDX.ID'].to_list()
            new_models = []
            for model in PDXNet_model_list:
                if model not in old_models:
                    new_models.append(model)
            if len(new_models) > 0:
                print("There are new models from PDXNet!!\n\n")
                print(new_models)
                self.flag = True



    def transform(self):
        return False
