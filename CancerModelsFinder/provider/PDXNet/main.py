from CancerModelsFinder.classes import API
from CancerModelsFinder.provider.PDXNet.resources import *

from os.path import join, exists
import pandas as pd



class PDXNet_API():
    def __init__(self):
        self.API = API()
    def request(self):
        response = self.API.request()


class PDXNet_data_transformations():
    file_paths = file_paths
    def __init__(self):
        self.ref_data = pd.read_csv(self.file_paths['JAX_endpoint_reference_file'], sep='\t')

    def process_models(self):
        model = pd.read_csv(self.file_paths['model']).drop('Unnamed: 0', axis=1)
        patient = pd.read_csv(self.file_paths['patient']).drop('Unnamed: 0', axis=1)
        sample = pd.read_csv(self.file_paths['sample']).drop('Unnamed: 0', axis=1)
        PDMR_mol = pd.read_csv(self.file_paths['PDMR_mol']).drop('Unnamed: 0', axis=1)
        PDTC_mol = pd.read_csv(self.file_paths['PDTC_mol']).drop('Unnamed: 0', axis=1)
        self.PDXNet_models = [model, patient, sample, PDMR_mol, PDTC_mol]

    def models_in_CM(self):
        models_in_cancerModels = pd.DataFrame()
        for provider in ['HCI-BCM', 'MDAnderson', 'WISTAR-MDAnderson-Penn', 'WUSTL']:
            models = pd.read_csv(join(self.file_paths['source'], provider, provider+'_metadata-pdx_model.tsv'), sep='\t')
            models = preprocess_templates(models)
            models['Provider'] = provider
            models_in_cancerModels = pd.concat([models_in_cancerModels, models])
        self.model_list_in_CM = models_in_cancerModels

    def compare_models(self):
        self.models_in_CM()
        old_models = self.model_list_in_CM['model_id'].to_list()
        PDXNet_model_list = self.PDXNet_models[1]['ContributorPDX.ID'].to_list()
        new_models = []
        for model in PDXNet_model_list:
            if model not in old_models:
                new_models.append(model)
        if len(new_models) > 0:
            print("There are new models from PDXNet!!\n\n")
            self.flag = True
        else:
            self.flag = False

    def transform(self):
        return False
    def main(self):
        self.process_models()
        self.compare_models()
        if self.flag:
            self.transform()


