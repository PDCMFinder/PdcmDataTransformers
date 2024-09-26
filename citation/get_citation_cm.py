import multiprocessing
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from datetime import datetime
from Bio import Entrez
from tqdm import tqdm

import requests
import xml.etree.ElementTree as ET

import requests
import xml.etree.ElementTree as ET
import re


# Function to fetch full text XML and find the sentence containing the citation
def fetch_and_find_citation_sentence(pmc_id, cited_pmid):
    # Construct the URL to fetch the full text XML
    pmc_id = get_pmcid_from_pmid(pmc_id)
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc_id}/fullTextXML"

    try:
        # Fetch the XML content from Europe PMC
        response = requests.get(url)
        response.raise_for_status()  # Check for request errors

        # Parse the XML content
        root = ET.fromstring(response.content)

        # Find all <xref> tags which usually contain citation information
        body_text = ""
        for body in root.findall('.//body'):
            body_text = ET.tostring(body, encoding='utf8', method='text').decode('utf8')

        # Create a pattern to search for the PMID (in citation format, e.g., [PMID:12345678])
        pattern = re.compile(rf"\[PMID:{cited_pmid}\]")

        # Split the body text into sentences
        sentences = re.split(r'(?<=[.!?]) +', body_text)

        # Search each sentence for the PMID
        for sentence in sentences:
            if pattern.search(sentence):
                return sentence.strip()  # Return the sentence containing the citation

        # If no sentence is found with the PMID
        return f"No sentence found with PMID {cited_pmid} in the citing article."

    except requests.exceptions.RequestException as e:
        return f"Error fetching full text for PMC ID {pmc_id}: {e}"



def get_pmcid_from_pmid(pmid):
    """
    Function to get PMCID from a given PubMed ID (PMID) using Europe PMC API.

    Args:
    pmid (str or int): The PubMed ID for which you want to retrieve the PMCID.

    Returns:
    pmcid (str): The corresponding PMCID, or None if no PMCID is found.
    """
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}&format=json"

    try:
        response = requests.get(url)
        data = response.json()

        # Check if the result exists and has a PMCID
        if data["hitCount"] > 0 and "pmcid" in data["resultList"]["result"][0]:
            return data["resultList"]["result"][0]["pmcid"]
        else:
            return None  # No PMCID found
    except Exception as e:
        print(f"Error fetching PMCID for PMID {pmid}: {e}")
        return None

# Function to fetch and process full text XML from Europe PMC and check for multiple model IDs
def fetch_and_search_models_in_full_text(pm_id, model_ids):
    # Construct the Europe PMC full-text XML URL
    pmc_id = get_pmcid_from_pmid(pm_id)
    if pmc_id is None:
        return False
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc_id}/fullTextXML"

    try:
        # Fetch the XML content from Europe PMC
        response = requests.get(url)
        response.raise_for_status()  # Raise error if the request fails

        # Parse the XML content
        root = ET.fromstring(response.content)

        # Find the <body> tag where the full text is located
        body_text = ""
        for body in root.findall('.//body'):
            # Extract all text within the body tag
            body_text = ET.tostring(body, encoding='utf8', method='text').decode('utf8')

        # Check if any model ID is mentioned in the full text
        for model_id in model_ids:
            if model_id in body_text:
                print(f"Model ID {model_id} found in full text of {pmc_id}.")
                return True

        print(f"None of the model IDs were found in full text of {pmc_id}.")
        return False

    except requests.exceptions.RequestException as e:
        print(f"Error fetching full text for PMC ID {pmc_id}: {e}")
        return False


def fetch_title_and_date(email, pmid):
    Entrez.email = email  # Always provide an email for NCBI API
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
        records = Entrez.read(handle)
        title = records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']

        # Extract publication date
        pub_date = records['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
        if 'Year' in pub_date and 'Month' in pub_date:
            pub_date_str = f"{pub_date['Year']}-{pub_date['Month']}-01"
        else:
            pub_date_str = f"{pub_date['Year']}-Jan-01"
        pub_date = datetime.strptime(pub_date_str, '%Y-%b-%d')
    except Exception as e:
        print(f"Error fetching data for PMID {pmid}: {e}")
        title, pub_date = None, None
    return title, pub_date


def fetch_citing_pmids(email, pmid):
    Entrez.email = email
    try:
        handle = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin")
        records = Entrez.read(handle)
        citing_pmids = [link['Id'] for link in records[0]['LinkSetDb'][0]['Link']]
    except Exception as e:
        print(f"Error fetching citing PMIDs for PMID {pmid}: {e}")
        citing_pmids = []
    return citing_pmids


class get_citations():
    def __init__(self):
        self.cpu_count = multiprocessing.cpu_count()
        self.email = "tushar@ebi.ac.uk"
        # Initialize lists to store results
        self.results = list()
        self.citations = list()
        # Set the filter date as April 1st, 2022
        self.filter_date = datetime(2022, 4, 1)

    def get_cm_metadata(self):
        model_publications = pd.read_json("https://www.cancermodels.org/api/model_information?publication_group_id=not.is.null")
        publication_id = pd.read_json("https://www.cancermodels.org/api/publication_group")
        model_publications = model_publications[['external_model_id', 'type', 'data_source', 'publication_group_id']].merge(publication_id, how='left', left_on='publication_group_id', right_on='id')
        df = model_publications[['external_model_id', 'data_source', 'pubmed_ids']].drop_duplicates()
        df = df.assign(pubmed_id=df['pubmed_ids'].str.replace(', ', ',').str.split(',')).explode('pubmed_id')
        df['pubmed_id'] = df['pubmed_id'].str.strip().str.replace('PMID: ', '').str.replace('PMID:', '')
        self.pmid_and_model_id = df.groupby('pubmed_id')['external_model_id'].agg(list).to_dict()
        df = model_publications[['data_source', 'pubmed_ids']].drop_duplicates()
        df = df.assign(pubmed_id=df['pubmed_ids'].str.replace(', ', ',').str.split(',')).explode('pubmed_id')
        df['pubmed_id'] = df['pubmed_id'].str.strip().str.replace('PMID: ', '').str.replace('PMID:', '')
        self.all_pmids = [l for l in list(set(df['pubmed_id'])) if l != '']

    def generate_result_dict(self, pmid):
        title, pub_date = fetch_title_and_date(self.email, pmid)
        citing_pmids = fetch_citing_pmids(self.email, pmid)

        # Filter out PMIDs that are already in the DataFrame
        filtered_citing_pmids = [citing_pmid for citing_pmid in citing_pmids if citing_pmid not in self.all_pmids]

        # Check publication date of each citing PMID
        model_in_cpmid = []
        citing_after_filter = []
        for citing_pmid in filtered_citing_pmids:
            if citing_pmid not in self.citations:
                self.citations.append(citing_pmid)
                _, citing_pub_date = fetch_title_and_date(self.email, citing_pmid)
                if citing_pub_date and citing_pub_date > self.filter_date:
                    citing_after_filter.append(citing_pmid)
                    model_in_cpmid.append({citing_pmid: fetch_and_find_citation_sentence(citing_pmid, pmid)})

        # Remove duplicates within the list
        citing_after_filter = list(set(citing_after_filter))

        # Append results
        self.results.append({
            'pubmed_id': pmid,
            'title': title,
            'cited_pmids': citing_after_filter,
            'citation_count': len(citing_after_filter),
            'cited_mention_text': model_in_cpmid,
        })

    def get_citations_for_cm(self):
        self.get_cm_metadata()
        for i in tqdm(self.all_pmids):
            self.generate_result_dict(i)
        #with ThreadPoolExecutor(max_workers=self.cpu_count) as executor:
        #    executor.map(self.generate_result_dict, self.all_pmids)
        if len(self.results) > 0:
            results_df = pd.DataFrame(self.results)
            codon_dir = "/hps/nobackup/tudor/pdcm/annotation-data/"
            results_df.to_csv(codon_dir+'citations_with_mentioned_lines.csv', index=False)
            #results_df['cited_pmids'] = results_df['cited_pmids'].apply(lambda x: list(set(x)))
            #all_cited_pmids = set(pm for sublist in results_df['cited_pmids'] for pm in sublist)
            #results_df['unique_citations'] = results_df['cited_pmids'].apply(lambda x: list(set(x) & all_cited_pmids))
            #results_df['unique_citation_count'] = results_df['unique_citations'].apply(len)
            #results_df.to_csv(codon_dir+'citations_with_ts.csv', index=False)


get_citations().get_citations_for_cm()