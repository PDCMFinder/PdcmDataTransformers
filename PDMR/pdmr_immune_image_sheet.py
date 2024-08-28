import os
import sys
from utils import read_metadata_without_fields, read_metadata_with_fields,  request
from os import getcwd
import pandas as pd
from tqdm import tqdm
import re
from bs4 import BeautifulSoup

def parse_html(response):
    soup = BeautifulSoup(response.text, 'html.parser')
    table = soup.find('table')
    data_rows = []
    for row in table.find_all('tr')[2:]:  # Skip the first two rows (headers and separator)
        data = [td.text.strip() for td in row.find_all('td')]
        if len(data)>0 and data[0].__contains__("General Model Information"):
            data_rows.append(data)
    m = re.search("MSI.*HLA.*\n.*\n.*", re.sub("\n+", "\n", data_rows[0][0]))
    if m:
        out = m.group()
        out = out.replace('MSI Status', 'MSI Status: ')
        out = out.replace("HLA Profile", "\nHLA Profile: ")
    else:
        out = ''
    if out.__contains__('MSI Status'):
        msi = out.split("\n")[0].split(" ")[2]
    else:
        msi = ''
    if out.__contains__("HLA Profile"):
        out = re.sub("\(s\)\n.*\n", "", out)
        out = re.sub("./.", "", out)
        hla = out.split('\n')[1].split(" ")[2]
    else:
        hla = ''
    return [msi, hla]

def get_immune_generate_files(home):
    ms_sheet = read_metadata_with_fields(os.path.join(home, "PDMR_molecular_metadata-sample.tsv"))
    immune_sheet = pd.DataFrame(
        columns=['Field', 'sample_id', 'marker_type', 'marker_name', 'marker_value', 'essential_or_additional_details',
                 'platform_id'])
    types = ["MMR", "TMB", "Mutations per mb", "Ploidy"]
    pid = {"MMR": 'immune_mmr', "TMB": 'immune_tmb', "Mutations per mb": 'immune_mpm', "Ploidy": 'immune_ploidy'}
    for i in tqdm(range(merged.shape[0]), "Fetch from PDMR"):
        mid = merged.iloc[i, 0]
        sid = merged.iloc[i, 1]
        so = merged.iloc[i, 2]
        link = merged.iloc[i, 3]
        rdu = merged.iloc[i, 4]
        response = request(link, False, 'get')
        dt = set()
        if response.text is not None:
            msi, hla = parse_html(response)
            if msi.__contains__("MSI"):
                msi = msi.replace("MSI-", "")
                ms_sheet = pd.concat([ms_sheet, pd.DataFrame([['', mid, sid, so, merged.iloc[i, 6], merged.iloc[i, 7],
                                                               merged.iloc[i, 8], merged.iloc[i, 9], rdu,
                                                               'immune_MSI']], columns=ms_sheet.columns)])
                row = [['', sid, 'Model Genomics', 'MSI', msi, '', 'immune_MSI']]
                immune_sheet = pd.concat([immune_sheet, pd.DataFrame(row, columns=immune_sheet.columns)])
                dt.add('MSI')
            if hla != '':
                pattern = re.compile(r'([A-Z]\*\d{2}:\d{2})')
                hla_types = pattern.findall(hla)
                hla_a_types = ['HLA-' + hla_type for hla_type in hla_types if hla_type.startswith('A')]
                hla_b_types = ['HLA-' + hla_type for hla_type in hla_types if hla_type.startswith('B')]
                hla_c_types = ['HLA-' + hla_type for hla_type in hla_types if hla_type.startswith('C')]
                hla_a_string = ', '.join(hla_a_types)
                hla_b_string = ', '.join(hla_b_types)
                hla_c_string = ', '.join(hla_c_types)
                row = [['', sid, 'HLA type', 'HLA-A', hla_a_string, '', 'immune_HLA_type'],
                       ['', sid, 'HLA type', 'HLA-B', hla_b_string, '', 'immune_HLA_type'],
                       ['', sid, 'HLA type', 'HLA-C', hla_c_string, '', 'immune_HLA_type']]
                immune_sheet = pd.concat([immune_sheet, pd.DataFrame(row, columns=immune_sheet.columns)])
                ms_sheet = pd.concat([ms_sheet, pd.DataFrame([['', mid, sid, so, merged.iloc[i, 6], merged.iloc[i, 7],
                                                               merged.iloc[i, 8], merged.iloc[i, 9], rdu,
                                                               'immune_HLA_type']], columns=ms_sheet.columns)])
                dt.add('HLA')
        for t in types:
            if t not in dt:
                dt.add(t)
                ms_sheet = pd.concat([ms_sheet, pd.DataFrame([['', mid, sid, so, merged.iloc[i, 6], merged.iloc[i, 7],
                                                               merged.iloc[i, 8], merged.iloc[i, 9], rdu, pid[t]]],
                                                             columns=ms_sheet.columns)])
                row = [['', sid, 'Model Genomics', t, 'Not provided', '', pid[t]]]
                immune_sheet = pd.concat([immune_sheet, pd.DataFrame(row, columns=immune_sheet.columns)])

    immune_sheet.to_csv(os.path.join(home, "PDMR_immunemarker-Sheet1_new.tsv"), sep='\t', index=False)
    ms_sheet.drop_duplicates()
    ms_sheet.to_csv(os.path.join(home, "PDMR_molecular_metadata-sample_2.tsv"), sep='\t', index=False)


def generate_image_sheet(model_metadata):
    images = pd.DataFrame()
    for i in tqdm(range(model_metadata.shape[0]), "Fetch from PDMR"):
        originator_df = pd.DataFrame()
        pdx_df = pd.DataFrame()
        invitro_df = pd.DataFrame()
        default_url = "https://pdmdb.cancer.gov/web/apex/"
        mid = model_metadata.iloc[i, 0]
        type = model_metadata.iloc[i, 1]
        link = model_metadata.iloc[i, 2]
        response = request(link, False, 'get')
        try:
            table = get_table_of_samples(response)
            ## originator
            originator_url = table[table['sampleID'] == 'ORIGINATOR'].reset_index(drop=True)['url']
            if originator_url.shape[0] > 0:
                originator_url = default_url + originator_url[0]
                originator_df = pd.concat(
                    [originator_df, process_originator_images(request(originator_url, False, 'get'), mid)])
            if type == "PDX":
                ## PDX
                pdx = table[table['pdm_type'] == 'PDX'].reset_index(drop=True)
                pdx['url'] = default_url + pdx['url']
                pdx_df = process_pdx_images(pdx, mid)
            if type != "PDX":
                sampleID = mid.split('-', 3)[-1]
                invitro = table[table['sampleID'] == sampleID].reset_index(drop=True)
                invitro['url'] = default_url + invitro['url']
                invitro_df = process_invitro_images(invitro, mid)
                invitro_df['type'] = type
                invitro_df['Field'] = ""
            temp = generate_image_df(originator_df, pdx_df, invitro_df)
            images = pd.concat([images, temp]).reset_index(drop=True)
        except Exception as e:
            # print(f"error for {mid}: {e}")
            continue
    return images


start_dir = getcwd()
if len(sys.argv) > 1:
    home = sys.argv[1]
    if os.path.exists(home):
        model_sharing = read_metadata_without_fields(os.path.join(home, "PDMR_metadata-sharing.tsv"))
        mol_sample = read_metadata_without_fields(os.path.join(home, "PDMR_molecular_metadata-sample.tsv"))
        merged = model_sharing.merge(mol_sample, on="model_id", how="left")
        merged = merged[merged.fillna(0)['sample_origin'] != 0]
        merged = merged[
            ["model_id", "sample_id", "sample_origin", "database_url", 'raw_data_url', "passage", 'host_strain_name',
             'host_strain_nomenclature', 'engrafted_tumor_collection_site', 'platform_id']]
        merged = merged.drop_duplicates(['model_id', 'sample_id']).reset_index(drop=True)
        get_immune_generate_files(home)
        model_metadata = pd.read_json(
            f"https://www.cancermodels.org/api/model_metadata?data_source=eq.PDMR&select=model_id,type,source_database_url")
        image_df = generate_image_sheet(model_metadata)
        image_df.to_csv(os.path.join(home, "PDMR_metadata-model_image.tsv"), sep='\t', index=False)


def get_table_of_samples(response):
    soup = BeautifulSoup(response.text, 'html.parser')
    data_rows = []
    for div in soup.find_all('div'):
        attrs = div.attrs
        if 'class' in attrs.keys() and attrs['class'][0] == 'rc-title':
            header = div.text.strip()
            if header == "Sample (PDX)":
                table = div.findNext().find_all('table')[0]
                for row in table.find_all('tr')[2:]:  # Skip the first two rows (headers and separator)
                    # print(row.find_all('td'))
                    data = [td.find('a').attrs['href'] if td.find('a') else td.text.strip() for td in row.find_all('td')
                            if 'class' in td.attrs.keys() and td.attrs['class'][0] == 'data']
                    if len(data) > 0:
                        data_rows.append(data)
    return pd.DataFrame(data_rows,
                        columns=['url', 'pdm_type', 'sampleID', 'patient_origin', 'pdx_passage', 'images', 'mutation',
                                 'wes', 'rnaseq'])


def process_originator_images(response, mid):
    soup = BeautifulSoup(response.text, 'html.parser')
    divs = [div.findNext().find_all('div') for div in soup.find_all('div') if
            'class' in div.attrs.keys() and div.attrs['class'][0] == 'rc-title' and div.text.strip() == 'Sample']
    table = [div.findNext().find_all('table') for div in divs[0] if
             'class' in div.attrs.keys() and div.attrs['class'][
                 0] == 'rc-title' and div.text.strip() == 'Pathology Data'][
        0][0]
    tt = pd.DataFrame([[td.find('a').attrs['href'] if td.find('a') else td.find('img').attrs['src'] if td.find(
        'img') else td.text.strip() for row in table.find_all('tr') for td in row.find_all('td') if
                        'class' in td.attrs.keys() and td.attrs['class'][0] == 'data'][0:8]],
                      columns=["View", "TumorGrade", "TumorContent", "Necrosis", "Stromal", "InflammatoryCells",
                               "Low_res_image", "High_res_image"])
    tt['sample_id'] = "ORIGIN"
    tt['model_id'] = mid
    tt['passage'] = "ORIGIN"
    return tt


def process_pdx_images(df, mid):
    temp = pd.DataFrame()
    for index, row in df.iterrows():
        response = request(row['url'], False, 'get')
        soup = BeautifulSoup(response.text, 'html.parser')
        divs = [div.findNext().find_all('div') for div in soup.find_all('div') if
                'class' in div.attrs.keys() and div.attrs['class'][0] == 'rc-title' and div.text.strip() == 'Sample']
        table = [div.findNext().find_all('table') for div in divs[0] if
                 'class' in div.attrs.keys() and div.attrs['class'][
                     0] == 'rc-title' and div.text.strip() == 'Pathology Data'][
            0][0]
        df_list = [td.find('a').attrs['href'] if td.find('a') else td.find('img').attrs['src'] if td.find(
            'img') else td.text.strip() for row in table.find_all('tr') for td in row.find_all('td') if
                   'class' in td.attrs.keys() and td.attrs['class'][0] == 'data']
        df_list = [df_list[i:i + 8] for i in range(0, len(df_list), 8)]
        tt = pd.DataFrame(df_list,
                          columns=["View", "TumorGrade", "TumorContent", "Necrosis", "Stromal", "InflammatoryCells",
                                   "Low_res_image", "High_res_image"])
        tt['sample_id'] = row['sampleID']
        tt['model_id'] = mid
        tt['passage'] = row['pdx_passage']
        temp = pd.concat([temp, tt])
    return temp.drop_duplicates()


def list_to_dataframe(flat_list):
    # Step 1: Convert flat list to 2D list with 4 elements per sub-list
    two_d_list = [flat_list[i:i + 4] for i in range(0, len(flat_list), 4)]

    # Step 2: Convert 2D list to DataFrame
    df = pd.DataFrame(two_d_list, columns=['Column1', 'Column2', 'Column3', 'Column4'])

    return df


def process_invitro_images(df, mid):
    temp = pd.DataFrame()
    for index, row in df.iterrows():
        response = request(row['url'], False, 'get')
        soup = BeautifulSoup(response.text, 'html.parser')
        divs = [div.findNext().find_all('div') for div in soup.find_all('div') if
                'class' in div.attrs.keys() and div.attrs['class'][0] == 'rc-title' and div.text.strip() == 'Sample']
        table = [div.findNext().find_all('table') for div in divs[0] if
                 'class' in div.attrs.keys() and div.attrs['class'][
                     0] == 'rc-title' and div.text.strip() == 'In Vitro Images'][
            0][0]
        df_list = [td.find('a').attrs['href'] if td.find('a') else td.find('img').attrs['src'] if td.find(
            'img') else td.text.strip() for row in table.find_all('tr') for td in row.find_all('td') if
                   'class' in td.attrs.keys() and td.attrs['class'][0] == 'data']
        df_list = [df_list[i:i + 4] for i in range(0, len(df_list), 4)]
        tt = pd.DataFrame(df_list, columns=["View", "image_type", "notes", "image"])
        tt['sample_id'] = row['sampleID']
        tt['model_id'] = mid
        tt['passage'] = row['pdx_passage']
        temp = pd.concat([temp, tt])
    return temp.drop_duplicates()


def melt_df(df):
    df_melted = df.melt(
        id_vars=['View', 'TumorGrade', 'TumorContent', 'Necrosis', 'Stromal', 'InflammatoryCells', 'sample_id',
                 'model_id', 'passage'], value_vars=['Low_res_image', 'High_res_image'], var_name='magnification',
        value_name='URL')
    df_melted['URL'] = "https://pdmdb.cancer.gov/web/apex/" + df_melted['URL']
    df_melted['magnification'] = df_melted['magnification'].str.replace('Low_res_image',
                                                                        'Low Magnification Image').str.replace(
        'High_res_image', 'High Magnification Image')
    df_melted['Field'] = ""
    return df_melted


def generate_image_df(og, pdx, invitro):
    if og.shape[0] > 0:
        og = melt_df(og)
        og['description'] = og.apply(lambda
                                         x: f"Tumor Grade: {x['TumorGrade']}, Tumor Content: {x['TumorContent']}, Necrosis: {x['Necrosis']}, Stromal: {x['Stromal']}, Inflammatory Cells: {x['InflammatoryCells']}, Sample Type: Patient, Staining: H&E",
                                     axis=1)
        og['sample_type'] = "patient"
        og['passage'] = ""
        og['staining'] = "H&E"
        og = og[['Field', 'model_id', 'URL', 'description', 'sample_type', 'passage', 'magnification', 'staining']]
    else:
        og = pd.DataFrame(
            columns=['Field', 'model_id', 'URL', 'description', 'sample_type', 'passage', 'magnification', 'staining'])
    if pdx.shape[0] != 0:
        pdx = melt_df(pdx)
        pdx['description'] = pdx.apply(lambda
                                           x: f"Tumor Grade: {x['TumorGrade']}, Tumor Content: {x['TumorContent']}, Necrosis: {x['Necrosis']}, Stromal: {x['Stromal']}, Inflammatory Cells: {x['InflammatoryCells']}, Sample Type: Patient, Staining: H&E",
                                       axis=1)
        pdx['sample_type'] = "xenograft"
        pdx['staining'] = "H&E"
        og = pd.concat([og, pdx[og.columns]]).reset_index(drop=True)
    if invitro.shape[0] != 0:
        invitro['URL'] = "https://pdmdb.cancer.gov/web/apex/" + invitro['image']
        invitro['description'] = invitro.apply(lambda x: f"Staining: {x['image_type']}, Pathology notes: {x['notes']}",
                                               axis=1)
        invitro['sample_type'] = invitro['type']
        invitro['staining'] = invitro['image_type']
        invitro['passage'] = ""
        invitro['magnification'] = ""
        og = pd.concat([og, invitro[og.columns]]).reset_index(drop=True)
    return og