class data_transformations():

    def preprocess_templates(self):
        if 'Field' in metadata.columns:
            metadata = metadata.loc[metadata.Field.astype('str').str.startswith('#') != True,].reset_index(drop=True)
            metadata = metadata.drop('Field', axis=1)
        return metadata