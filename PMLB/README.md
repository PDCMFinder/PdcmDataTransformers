#Retrieving PMLB omic data

1. You need an account on cbioportal.ca and sufficient permissions to access 
the omic data. Be sure that the contact at PMLB as given you permissions.

2. Log in to cbioportal.ca and then go to https://cbioportal.ca/cbioportal/api/swagger-ui.html by entering the url into your browser. The swagger-ui is in beta and is has no obvious link from the cbioportal.ca.

3. 
    - Click /molecular-profiles
    - Click 'try it out' and then 'Execute'
The output will show all molecular profiles that are available for your account. If your desired profile is not present for you then you likely do not have permissions. Contact PMLB. Get both the molecular profile and study_id

4. 
    - Click on either the Mutations or the Discrete Copy Number Alterations button on the swagger api. 
    - Then click the Get request.
    - Click 'try it out'
    - Either data types will require the molecular profile id and the Sample List ID. The sample list ID is the study_id with "_all" appended to it. 
    - Set the detail of the response to "detailed     - In CNA the copy number event can be set to 'ALL'

5. 
    - Hit execute. If the data set is small you can hit 'download' else use the generated Request URL if the data set is large. If the data is large it can take more than a few minutes to load it into the website. Click download. This provides the JSON neccessary to then convert to the PDCM Finder mutation or cna format  
