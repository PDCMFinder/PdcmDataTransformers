import pandas as pd
class data_transformations():
    def __init__(self):
        self.metadata = pd.DataFrame()

    def drop_template_rows_columns(self, metadata = pd.DataFrame()):
        if 'Field' in metadata.columns:
            metadata = metadata.loc[metadata.Field.astype('str').str.startswith('#') != True,].reset_index(drop=True)
            metadata = metadata.drop('Field', axis=1)
        return metadata