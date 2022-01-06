    "The dwarves delved too greedily and too deep..."
    
This project transforms the HCMI data set into the PDCM Finder format. There are two main proccesses:
    1. Parse the metadata
    2. Download and parse somatic mutations.

To parse metadata two files are required:
    1. model-table.tsv This can be downloaded from https://hcmi-searchable-catalog.nci.nih.gov/ Click export on the top right of the table and click "TSV(all columns)"
    2. clinical.tsv is in a .tar.gz file that can be found at the GDC repository for the project. https://portal.gdc.cancer.gov/projects/HCMI-CMDC Click the clinical button at the top right of the page. Then "tar -xvf" the file to extract the content. Place the clinical.tsv into the resource directory.

Templates are required to run the script, currently. However this is soon to be deprecated. The exist in the resource directory.

To run HcmiToPdcm.py you will need python3 and pandas installed.

The script will output:
    - HCMI_metadata-patient.tsv
    - HCMI_metadata-patient_sample.tsv
    - HCMI_metadata-cell_model.tsv
    - HCMI_metadata-sharing.tsv

Because of the absence of validation data no tsv is created. It will need to be manually created.

Next to download and parse the mutation data run the script: get_somatic_files.py 

This uses the id's found in the model-table.tsv to download maf.gz's from the GDC API. The are unzipped and transformed to the PDCM format. Both the HCMI_molecular_metadata-sample.tsv and the HCMI_molecular_metadata-platform.tsv. 


