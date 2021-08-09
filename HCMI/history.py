import pandas as pd
import sys
import os
os.chdir("/home/afollette/Jupyter_notebooks/HCMI/")
hcmiDf = pd.read_excel("model-table.tsv", engine="openpyxl")
hcmiDf = pd.read_tsv("model-table.tsv", sep='\t')
hcmiDf = pd.read_csv("model-table.tsv", sep='\t')
atient_template = pd.read_csv("metadata_template-patient.tsv", sep='\t').dropna(axis=0, how='all')
sample_template = pd.read_csv("metadata_template-sample.tsv", sep='\t').dropna(axis=0, how='all')
%history -f history.py
