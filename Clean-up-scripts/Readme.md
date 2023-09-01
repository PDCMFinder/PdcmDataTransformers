These scripts are for cleaning up the data in the CancerModels.org repository.

- Data_cleanup JupyterNotebook:
  Can perform the following tasks
  - Generates XLSX files for each provider
  - Fixes typos in collection site
  - Convert data type from string to number for molecular data
  - Generates File list for submission to BIA
- Generate data cBioPortal JupyterNotebook:
  Generates the clinical and molecular data files for cBioPortal instance. 
- Merge All data to One JupyterNotebook:
  Generates one Excel file will all the data
- Project to PROVIDER JupyterNotebook:
  Splits individual projects such as PIVOT, HCMI to individual provider based on the source of the model
- Remove empty Columns and rows pythons script:
  Removes empty columns and rows from each sheet in the data repo.
