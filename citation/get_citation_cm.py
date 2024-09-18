from tqdm import tqdm
import pandas as pd
from datetime import datetime
from Bio import Entrez


def fetch_title_and_date(pmid):
    Entrez.email = "tushar@ebi.ac.uk"  # Always provide an email for NCBI API
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


# Function to fetch citing publications by PMID
def fetch_citing_pmids(pmid):
    Entrez.email = "tushar@ebi.ac.uk"
    try:
        handle = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin")
        records = Entrez.read(handle)
        citing_pmids = [link['Id'] for link in records[0]['LinkSetDb'][0]['Link']]
    except Exception as e:
        print(f"Error fetching citing PMIDs for PMID {pmid}: {e}")
        citing_pmids = []
    return citing_pmids

model_publications = pd.read_json("https://www.cancermodels.org/api/model_information?publication_group_id=not.is.null")
publication_id = pd.read_json("https://www.cancermodels.org/api/publication_group")
model_publications = model_publications[['external_model_id', 'type', 'data_source', 'publication_group_id']].merge(publication_id, how='left', left_on='publication_group_id', right_on='id')
df = model_publications[['data_source', 'pubmed_ids']].drop_duplicates()
df = df.assign(pubmed_id=df['pubmed_ids'].str.replace(', ', ',').str.split(',')).explode('pubmed_id')
df['pubmed_id'] = df['pubmed_id'].str.strip().str.replace('PMID: ', '').str.replace('PMID:', '')
df.to_csv('data_source_publications.csv', index=False)
# Collect all PMIDs in the data
all_pmids = set(df['pubmed_id'])

# Initialize lists to store results
results = []
citations = list()
# Set the filter date as April 1st, 2022
filter_date = datetime(2022, 4, 1)

# Process each PMID
for i in tqdm(range(1, len(all_pmids))):
    pmid = list(all_pmids)[i]
    title, pub_date = fetch_title_and_date(pmid)
    citing_pmids = fetch_citing_pmids(pmid)

    # Filter out PMIDs that are already in the DataFrame
    filtered_citing_pmids = [citing_pmid for citing_pmid in citing_pmids if citing_pmid not in all_pmids]

    # Check publication date of each citing PMID
    citing_after_filter = []
    for citing_pmid in filtered_citing_pmids:
        if citing_pmid not in citations:
            citations.append(citing_pmid)
            _, citing_pub_date = fetch_title_and_date(citing_pmid)
            if citing_pub_date and citing_pub_date > filter_date:
                citing_after_filter.append(citing_pmid)

    # Remove duplicates within the list
    citing_after_filter = list(set(citing_after_filter))

    # Append results
    results.append({
        'pubmed_id': pmid,
        'title': title,
        'cited_pmids': citing_after_filter,
        'citation_count': len(citing_after_filter)
    })


# Create DataFrame from results
results_df = pd.DataFrame(results)

# Remove duplicates within 'cited_pmids'
results_df['cited_pmids'] = results_df['cited_pmids'].apply(lambda x: list(set(x)))

# Create a set of all unique cited PMIDs across the DataFrame
all_cited_pmids = set(pm for sublist in results_df['cited_pmids'] for pm in sublist)

# Create 'unique_citations' column with unique cited PMIDs
results_df['unique_citations'] = results_df['cited_pmids'].apply(lambda x: list(set(x) & all_cited_pmids))

# Add a new column for the count of unique citations
results_df['unique_citation_count'] = results_df['unique_citations'].apply(len)

results_df.to_csv('citations_updated.csv', index=False)