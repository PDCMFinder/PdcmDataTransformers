This repository is used to transform the Cell Model Passports data into the PDCM Finder
format.

The contact points for data at Cell Model Passports was Syd contacted.

1. The script requires you retrieve the data from the Cell Model Passports website. Go here https://cellmodelpassports.sanger.ac.uk/downloads and click to download the annotations for the modesl. Previously the file was called model_list_20210310.xlsx.
2. Convert and save file to csv.
3. Rename variable in cellPassportToPDCMFinder.py to read the model annotations.
4. To run the cellPassportToPDCMFinder.py script use python3 with pandas installed.
5. Templates are required to run the pipeline and are supplied in this directory. Changes these templates as neccessary.

6. Files are ouput into the pdcm_format directory. These are the medata sheets:
    * patient
    * patient sample
    * cell model
    * validation
    * sharing

Post-proccessing steps:
    
    1. Validation data was recieved in a separate file. This file was commited to this repo as "AuthenticatedModelList_Organoids-model_list_20210719.tsv". This information should be manually added to the metadata-validation.tsv sheet and quality controlled. Potentially it contains models that did not make it to the transformed metadata. These should be removed. As per communication the validaiton technology should be "SNP"

    2. Web content was already prepared and is in the web directory.  

    3. Molecular data exists, but has yet to be loaded. If this is to be done the Moleculare_metadata templates will need to be created.

