import pandas as pd
import numpy as np
import re
import sys
import string
import ipdb
from fuzzywuzzy import process
import argparse

#This file adds new metadata to the NCBI Virus metadata
#The files get checked and curated if necessary

# Function to check if an accession number is real: it uses the entrez functionality of ncbi
from check_accession import extract_accession

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add additional metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  metavar=' ', help="input metadata")
    parser.add_argument('-o', '--output', metavar=' ', help="output metadata")
    parser.add_argument('--add', help="input additional data file")
    parser.add_argument('--local', help="input local accession file")
    parser.add_argument('--regions', help="file to specify regions: format = country region")
    parser.add_argument('--id', help="id: strain or accession", choices=["strain","accession"],default="accession")
    parser.add_argument('--rename', help="copy of updated_metadata.tsv")
    parser.add_argument('--C1like', help="accessions to C1-like subgenotypes")
    parser.add_argument('--update', help="date when sequences were added")
    args = parser.parse_args()

    id_field = args.id
    input_csv_meta = args.input
    output_csv_meta = args.output
    add_data = args.add # if several files, use more than one assignment (add_data_1, add_data_2,...)
    con_reg_table = args.regions
    local_accn = args.local
    renamed_strains = args.rename
    C1like_accn = args.C1like
    last_updated_file = args.update

    # load data
    meta = pd.read_csv(input_csv_meta, keep_default_na=True, sep='\t', index_col=False)
    new_data = pd.read_csv(add_data, keep_default_na=True, sep='\t', index_col=False)
    local_accn_file= pd.read_csv(local_accn, keep_default_na=True, sep='\t', index_col=False)
    renamed_strains_df = pd.read_csv(renamed_strains, keep_default_na=True, sep='\t', index_col=False,names=["accession","strain"])
    last_updated=pd.read_csv(last_updated_file, keep_default_na=True, sep='\t', index_col=False,names=["accession","date_added"])
    
    # # Identify records that need updating
    # needs_update = meta[~meta['accession'].isin(renamed_strains_df['accession'])]
    # needs_update.to_csv("data/no_strain_correction.tsv", sep='\t', index=False)

    # Create a lookup dictionary for strain updates
    lookup_strain = renamed_strains_df.set_index('accession')['strain'].to_dict()

    # Update strains in metadata according to lookup
    meta['strain'] = meta.apply(
        lambda row: lookup_strain.get(row['accession'], row['strain']),
        axis=1
    )

    # Create a dictionary for quick lookup
    accession_dict = local_accn_file.set_index('sample_name')['seq_accession'].to_dict()

    # Replace missing accessions using the dictionary
    new_data['accession'] = new_data.apply(
        lambda row: accession_dict.get(row['strain'], row['accession']),
        axis=1
    )

    # replace None values with NaN
    new_data['accession'] = new_data['accession'].replace([None], np.nan)

                
    ## Remove duplicates based on id_field
    new_data = new_data.drop_duplicates(subset=id_field)
    meta = meta.drop_duplicates(subset=id_field)

    # If location is missing, replace it with division
    meta['location'] = meta['location'].mask(meta['location'].isna(), meta['division'])

    # step 1: merge both files with to accession number
    new_meta = pd.merge(meta, new_data, on=id_field, how='outer').dropna(subset="accession")

    # add date_added column
    new_meta= pd.merge(new_meta,last_updated, on=id_field,how='left')

    # Creating the new strain column based on the conditions
    new_meta['strain'] = new_meta['strain_y'].mask(new_meta['strain_y'] == new_meta['accession'], new_meta['strain_x'])  # Take strain_x if strain_y == accession
    new_meta['strain'] = new_meta['strain'].mask(new_meta['strain_x'] == new_meta['accession'], new_meta['strain_y'])  # Take strain_y if strain_x == accession
    new_meta['strain'] = new_meta['strain'].mask(new_meta['strain_x'].isna(), new_meta['strain_y'])  # Take strain_y if strain_x is NaN
    new_meta['strain'] = new_meta['strain'].fillna(new_meta['strain_x'])  # Take strain_x if strain_y is NaN

    # Keep only the dates from assign_publications.tsv table - except if they're NA
    new_meta['date_new'] = new_meta['collection_date'].mask(sum([(new_meta['collection_date'].isna()),(new_meta['collection_date'] == "XXXX-XX-XX")])>=1, new_meta['date'])  # If date_y is unknown, keep date_x

    # Keep date_y if date_correction specified in origin
    new_meta['date_new'] = new_meta['date_new'].mask(new_meta['origin']=='date_correction', new_meta['collection_date'])

    new_meta['date'] = new_meta['date_new']

    # Region: keep to most detailed one (longest string)
    new_meta['region'] = new_meta['region_x'].mask(new_meta['region_x'].isna(), new_meta['region_y'])
    new_meta['region'] = new_meta['region'].mask(new_meta['region_x'].str.len()<new_meta['region_y'].str.len(), new_meta['region_y'])

    # Country: keep the non-missing ones
    new_meta['country'] = new_meta['country_y'].mask(new_meta['country_y'].isna(), new_meta['country_x'])
    new_meta['country'] = new_meta['country'].replace("Viet Nam", "Vietnam")
    new_meta['country'] = new_meta['country'].replace("UK", "United Kingdom")

    # Country: keep the non-missing ones
    new_meta['place'] = new_meta['place'].mask(new_meta['place'].isna(), new_meta['location'])

    # Clades: keep non-missing clades - subgenogroup
    new_meta['subgenogroup'] = new_meta['subgenogroup'].mask(new_meta['subgenogroup'].isna(), new_meta['clade_x'])

    # Isolation source: standardize
    # Function to map non-standard terms to standard terms
    def standardize_isolation_source(value):
        # Add mappings for non-standard terms to standard terms
        mapping = {
            'stool': 'feces',
            'Stool': 'feces',
            'Faeces': 'feces',
            'CSF': 'cerebrospinal fluid',
            'Cerebrospinal Fluid': 'cerebrospinal fluid',
            'Cerebrospinal Fluid ': 'cerebrospinal fluid',
            'Throat Swab': 'oronasopharynx',
            'Throat swab': 'oronasopharynx',
            'Oral Swab': 'oronasopharynx',
            'Mouth Swab': 'oronasopharynx',
            'Nasopharyngeal Aspirate': 'oronasopharynx',
            'Brain': 'brain',
            'Cerebellum': 'brain',
            'Left brain': 'brain',
            'Brainstem': 'brain',
            'Vesicle Swab': 'swab',
            'Ulcer Swab': 'lesion',
            'Rectal swab': 'swab',
            'Rectal Swab': 'swab',
            'Mouth': 'oronasopharynx',
            'Throat': 'oronasopharynx',
            'Tonsil': 'oronasopharynx',
            'Vesicle': 'lesion',
            'Spinal cord ': 'spinal cord'
        }
        
        val=mapping.get(value, value)
        if pd.isna(val):
            return val
        val=val.title()

        # Return the mapped value if it exists, otherwise return the original value
        return val
    
    # Apply the standardization to both columns
    new_meta['isolate-lineage-source'] = new_meta['isolate-lineage-source'].apply(standardize_isolation_source)
    new_meta['isolation'] = new_meta['isolation'].apply(standardize_isolation_source)

    # Combine the two columns, prioritizing the standardized values
    new_meta['combined_isolation_source'] = new_meta['isolation'].mask(new_meta['isolation'].isna(), new_meta['isolate-lineage-source'])

    # Replace the original isolation_source column
    new_meta['isolation_source'] = new_meta['combined_isolation_source']

    # Drop unnecessary columns
    new_meta = new_meta.drop(['isolation', 'combined_isolation_source'], axis=1)


    # Define a dictionary for old naming formats and spelling mistakes
    corrections = {
        'czech republic': 'czechia',
        'hongkong': 'hong_kong',
        'viet nam': 'vietnam',
        'uk': 'united kingdom',
        'ivory coast': 'côte d\'ivoire',
    }

    # Function to correct country names
    def correct_country_name(country):
        if pd.notna(country):
            country = country.strip().lower()
            corrected_country = corrections.get(country, country).title().replace('_', ' ')
            # print(f"Correcting '{country}' to '{corrected_country}'")  # Debugging statement
            return corrected_country
        return country

    # Apply the corrections to the 'country' column
    new_meta['country'] = new_meta['country'].apply(correct_country_name)

    # Check if the regions file is supplied
    if con_reg_table:
        # Read the regions file and create a dictionary for country-region mappings
        with open(con_reg_table) as f:
            regions = {line.split("\t")[0].strip().lower(): line.split("\t")[1].strip() for line in f.readlines()[1:]}

        # Function to get region from country
        def get_region(coun):
            if pd.notna(coun):
                coun = coun.strip().lower()
                return regions.get(coun, regions.get(coun.replace(' ', '_'), "NA")).replace('_', ' ').title()
            return "NA"

        # Update the 'region' column in the new_meta DataFrame with the new region values
        new_meta['region'] = new_meta['country'].apply(get_region)

    # Debugging statement to check the unique values in the 'country' column after correction
    print("Unique countries after correction:", new_meta['country'].unique())

    new_meta['has_diagnosis'] =~new_meta['diagnosis'].isna()

    # Define a mapping for full terms to their abbreviations and standardized names
    short_versions = {
        'poliomyelitis-like disease': 'AFP', 'acute flaccid paralysis': 'AFP',
        'paralysis': 'AFP', 'poliomyelitis-like paralysis': 'AFP', 'afp': 'AFP',
        'acute flaccid myelitis': 'AFP', 'afm': 'AFP', 'polio': 'AFP',
        'difficulty walking': 'AFP','paralyses': 'AFP',
        'Cerebellar ataxia': 'Ataxia',
        'Cerebellitis': 'Ataxia',

        'upper respiratory tract infection': 'Respiratory',
        'acute respiratory disease': "Respiratory",
        'asthmatic bronchitis': "Respiratory",
        'Cough': "Respiratory",
        'Influenza-like illness': "Respiratory",
        'Respiratory Symptoms':'Respiratory',

        'guillain-barré syndrome': 'GBS', 'guillain-barrésyndrome': 'GBS',

        'febrile illness': 'Fever', 'pyrexia': 'Fever', 'fever': 'Fever',

        'v&d': 'Vomiting; Diarrhea',

        'transverse myelitis': 'Myelitis',

        'diarrhoea': 'Diarrhea', 'diarrhea': 'Diarrhea',
        'Gastroenteritis':'GI infection','Gastrointestinal':'GI infection',

        'severe': "Severe HFMD", 
        'light': "Light HFMD", 'mild': "Mild HFMD",
        'Atypical hand foot and mouth disease': "Atypical HFMD",
        'atypical hfmd': "Atypical HFMD",
        'HFMD; Onychomadesis': "Atypical HFMD",
        'HFMD, Onychomadesis': "Atypical HFMD",
        'hfmd': 'HFMD', 'hand foot and mouth disease': 'HFMD',
        'hand-foot and mouth disease': 'HFMD',
        'Toe lesion':'Lesion',

        'Gingivostomatitis ':'Stomatitis',

        'sepsis': 'Sepsis','neonatal sepsis':'Sepsis',
        'Seizure':'Seizure',

        'herpangina': 'Herpangina', 'herpangine': 'Herpangina',

        'hypotension': 'Hypotension', 'hyperglycemia': 'Hyperglycemia',

        'severe pulmonary hemorrhage': 'Cardiopulmonary Collapse',
        'pulmonary hemorrhage': 'Cardiopulmonary Collapse',
        'Acute cardiogenic shock': 'Cardiogenic shock',
        'Acute cardiogenic shock, fatal': 'Cardiogenic shock; Fatality',
        'cardiopulmonary collapse': 'Cardiopulmonary Collapse',

        'vesicular': 'vesicular rash', 'rash with vesicles': 'vesicular rash',

        '(death)': 'Fatality', 'fatal': 'Fatality', 'death': 'Fatality', 'Died':'Fatality',
        'Shock':'Shock',

        'oral ulcers': 'Oral Ulcers',
        'Onychomadesis':'Onychomadesis',

        'cns symptoms': 'CNS', 'cns involvement': 'CNS', 'cns disorder': 'CNS',
        'neurological': 'CNS', 'severe neurological involvement': 'CNS',

        'myoclonic jerk': 'Myoclonic jerk', 'myoclonus': 'Myoclonic jerk',

        'lethargy': 'Lethargy', 'lethargic': 'Lethargy',

        'neck stiffness; vomiting': 'Meningitis',
        'aseptic meningitis': 'Aseptic Meningitis',
        'Benign intracranial hypertension':'Meningitis',

        'meningoencephalitis': 'Encephalitis', 'encephalytis': 'Encephalitis',
        'Encephalitis':'Encephalitis',
        'Encéphalite':'Encephalitis',

        'skin rash': 'Rash', 'exanthema': 'Rash', 'rash': 'Rash','Skin':'Rash',
        'eczema': 'Rash',
    }

    short_forms = set(short_versions.values())

    major_versions = {
        # CNS-related conditions
        'Encephalitis': 'CNS','Aseptic Meningitis': 'CNS','Meningitis': 'CNS','Opsomyoclonus Syndrome': 'CNS','Myoclonic jerk':'CNS',
        'CNS': 'CNS', 'lethargy':'CNS','Seizure':'CNS',

         # Fatality-related conditions
        'Fatality':'Fatality',

        # AFP-related conditions
        'GBS': 'AFP','Myelitis': 'AFP', 'AFP':'AFP', 'ataxia':'AFP',  

        # HFMD-related conditions
        'HFMD':'HFMD','Atypical HFMD':"Atypical HFMD", 'Onychomadesis':"Atypical HFMD",
        'Severe HFMD': 'HFMD', 'Light HFMD': 'HFMD', 'Mild HFMD': 'HFMD',

        # Other categories
        'Sepsis': 'other',
        'URTI': 'other',
        'Respiratory Symptoms': 'other',
        'Fever': 'other',
        'Vomiting; Diarrhea': 'other',
        'Diarrhea': 'other',
        'Oral Ulcers': 'other',
        'Hypotension': 'other',
        'Hyperglycemia': 'other',
        'Cardiogenic shock': 'other',
        'Cardiopulmonary Collapse': 'other',
        'Rash': 'other',
        'Vesicular Rash': 'other'   
    }

    major_forms = set(major_versions.values())

    def clean_diagnosis(diagnosis, threshold=75):
        if pd.isna(diagnosis):
            return np.nan
        
        # Check if the diagnosis is already a short form
        if diagnosis in short_forms:
            return diagnosis
        
        # Remove punctuation and split multiple diagnoses
        clean_diag = diagnosis.replace(',', ';').replace(' or ', ';').replace('/', ';').replace('  ', ' ').strip()
        diagnoses = [diag.strip() for diag in clean_diag.split(';')]

        # Standardize diagnoses and replace full terms with abbreviations
        standardized_diagnoses = []
        for diag in diagnoses:
            diag_lower = diag.lower()
            if diag_lower in short_versions:
                standardized_diagnoses.append(short_versions[diag_lower])
            else:
                # Use fuzzy matching to handle typos
                match = process.extractOne(diag_lower, short_versions.keys(), score_cutoff=threshold)
                if match:
                    standardized_diagnoses.append(short_versions[match[0]])
                else:
                    # Check if the original diagnosis is in short_forms
                    standardized_diagnoses.append(diag if diag in short_forms else diag.title())
        
        # Join the cleaned and standardized diagnoses back into a string
        return '; '.join(sorted(set(standardized_diagnoses)))

    def extract_major_diagnosis(cleaned_diagnosis, threshold=80):
        if pd.isna(cleaned_diagnosis) or cleaned_diagnosis == "":
            return np.nan
        
        # Check if the diagnosis is already a major category
        if cleaned_diagnosis in major_forms:
            return cleaned_diagnosis
        
        # Remove punctuation and split multiple diagnoses
        diagnoses = cleaned_diagnosis.split('; ')

        # Standardize diagnoses and replace full terms with major categories
        major_diagnoses = []
        for diag in diagnoses:
            if diag in major_versions:
                major_diagnoses.append(major_versions[diag])
        
        # Join the cleaned and standardized diagnoses back into a string
        return '; '.join(sorted(set(major_diagnoses)))

    # All the diagnosis
    new_meta['med_diagnosis_all'] = new_meta['diagnosis'].apply(lambda x: clean_diagnosis(x))

    # only the major diagnosis
    new_meta['med_diagnosis_major'] = new_meta['med_diagnosis_all'].apply(lambda x: extract_major_diagnosis(x))

    # Add filter for age and add age ranges
    new_meta['age_yrs'] = pd.to_numeric(new_meta['age_yrs'], errors='coerce')
    new_meta["has_age"] = ~new_meta["age_yrs"].isna()

    # parse gender
    new_meta['gender'] = new_meta['gender'].mask(new_meta['gender'].str.contains('female', case=False, na=False), 'F')
    new_meta['gender'] = new_meta['gender'].mask(new_meta['gender'].str.contains('male', case=False, na=False), 'M')

    #Define age bins and labels for years
    bins_years = [-np.inf, 1, 6, 11, 18, np.inf]
    labels_years = ['<=1 y/o', '1-5 y/o', '6-10 y/o', '11-17 y/o', '18+ y/o']

    # Define age bins and labels for months
    bins_months = [-np.inf, 0.25, 0.5, 1]
    labels_months = ['0-3 m/o', '4-6 m/o', '7-12 m/o']

    # Create age_range column using pd.cut for years
    new_meta['age_range'] = pd.cut(new_meta['age_yrs'], bins=bins_years, labels=labels_years, right=False).astype(str)

    # Handle ages less than 1 year old separately using pd.cut for months
    mask_months = new_meta['age_yrs'] < 1
    new_meta.loc[mask_months, 'age_range'] = pd.cut(new_meta.loc[mask_months, 'age_yrs'], bins=bins_months, labels=labels_months, right=False).astype(str)

    # rename length to NCBI_length_genome
    new_meta.rename(columns={"length": "NCBI_length_genome"}, inplace=True)

    # if ENPEN in origin, set ENPEN to True
    new_meta = new_meta.assign(ENPEN=new_meta['origin'].str.contains('ENPEN', case=False, na=False))

    # if ENPEN=TRUE; Authors in doi should be moved to 'authors'
    # e.g. Private: .... remove Private from authors, keep private in doi
    new_meta['authors'] = new_meta.apply(
        lambda row: row['doi'].replace("Private: ", "") if row['ENPEN'] else row['authors'],
        axis=1
    )
    new_meta['doi'] = new_meta.apply(
        lambda row: "Private" if row['ENPEN'] else row['doi'],
        axis=1
    )

    # write new metadata file to output
    new_meta= new_meta.loc[:,['accession', 'accession_version', 'strain', 'date', 'region', 'place',
        'country', 'host', 'gender', 'age_yrs','age_range',"has_age", 'has_diagnosis','med_diagnosis_all','med_diagnosis_major',
        'isolation_source', 'NCBI_length_genome',
        'subgenogroup','date_released',
         'abbr_authors', 'authors', 'institution','ENPEN','doi',
        'qc.overallScore', 'qc.overallStatus',
        'alignmentScore', 'alignmentStart', 'alignmentEnd', 'genome_coverage','date_added']]


    accn = pd.read_csv(C1like_accn, sep = "\t")

    # merge the two dataframes. Replace the subgenotype in meta with the subgenotype in accn
    new_meta2 = pd.merge(new_meta, accn, on="accession", how="outer")

    # Replace the subgenogroup and doi in meta with the subgenogroup and doi in accn if they are not NA
    new_meta2['subgenogroup'] = new_meta2['subgenogroup_y'].combine_first(new_meta2['subgenogroup_x'])

    # Replace the doi in meta with the doi in accn if they are not part of ENPEN
    new_meta2["doi"] = new_meta2["doi_x"].mask(new_meta2["ENPEN"] !=True, new_meta2["doi_y"])

    # Drop the temporary columns
    new_meta2.drop(columns=['subgenogroup_x', 'subgenogroup_y', 'doi_x', 'doi_y'], inplace=True)
    new_meta2 = new_meta2.drop_duplicates(subset="accession",keep="first")


    new_meta2.to_csv(output_csv_meta, sep='\t', index=False)