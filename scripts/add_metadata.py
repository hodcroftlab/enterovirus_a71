import pandas as pd
import numpy as np
import re
import sys
import string
import ipdb
from fuzzywuzzy import process
import argparse


#Is purely to add metadata when something new comes in last-minute without re-running everything
#Format of the file:
#column 1: 'strain' or 'accession' to indicate how to indicate the record to be changed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add additional metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  metavar=' ', help="input metadata")
    parser.add_argument('-o', '--output', metavar=' ', help="output metadata")
    parser.add_argument('--add', help="input additional data file")
    parser.add_argument('--regions', help="file to specify regions: format = country region")
    parser.add_argument('--id', help="id: strain or accession", choices=["strain","accession"],default="accession")
    args = parser.parse_args()

    id_field = args.id
    input_csv_meta = args.input
    output_csv_meta = args.output
    add_data = args.add # if several files, simply use more than one assignment (add_data_1, add_data_2,...)
    con_reg_table = args.regions

    # load data
    meta = pd.read_csv(input_csv_meta, keep_default_na=True, sep='\t', index_col=False)
    new_data = pd.read_csv(add_data, keep_default_na=True, sep='\t', index_col=False)
             
    ## Remove duplicates based on id_field
    # new_data = new_data.drop_duplicates(subset=id_field)
    # meta = meta.drop_duplicates(subset=id_field)

    # If location is missing, replace it with division
    meta['place'] = meta['location'].mask(meta['location'].isna(), meta['division'])

    # step 1: merge both files with to accession number
    new_meta = pd.merge(meta, new_data, on=id_field, how='left')

    # rename dupicated columns
    new_meta.rename(columns={"date":"date_x","collection_date":"date_y"},inplace=True)

    # Creating the new strain column based on the conditions
    new_meta['strain'] = new_meta['strain_y'].mask(new_meta['strain_y'] == new_meta['accession'], new_meta['strain_x'])  # Take strain_x if strain_y == accession
    new_meta['strain'] = new_meta['strain'].mask(new_meta['strain_x'] == new_meta['accession'], new_meta['strain_y'])  # Take strain_y if strain_x == accession
    new_meta['strain'] = new_meta['strain'].mask(new_meta['strain_x'].isna(), new_meta['strain_y'])  # Take strain_y if strain_x is NaN
    new_meta['strain'] = new_meta['strain'].fillna(new_meta['strain_x'])  # Take strain_x if strain_y is NaN

    # Remove "00:00:00" from the date strings
    new_meta['date_y'] = new_meta['date_y'].str.replace(' 00:00:00', '', regex=False)
        
    # Pick the row with the most accurate date
    new_meta['date'] = new_meta['date_x'].mask(new_meta['date_x'] == "XXXX-XX-XX", new_meta['date_y'])  # If date is unknown, replace with known
    new_meta['date'] = new_meta['date'].mask(new_meta['date_y'].isna(), new_meta['date_x'])  # If date_y is unknown, keep date_x

    # Function to count the number of 'XX' components in a date string
    def count_unknowns(date_str):
        return date_str.count('XX') if pd.notna(date_str) else float('inf')

    ## Keep the one with less XX
    new_meta['date'] = new_meta.apply(
        lambda row: row['date_y'] if (pd.notna(row['date_y']) and count_unknowns(row['date_y']) < count_unknowns(row['date_x'])) else row['date_x'], 
        axis=1
    )

    # Region: keep to most detailed one (longest string)
    new_meta['region'] = new_meta['region_x'].mask(new_meta['region_x'].isna(), new_meta['region_y'])
    new_meta['region'] = new_meta['region'].mask(new_meta['region_x'].str.len()<new_meta['region_y'].str.len(), new_meta['region_y'])


    # Country: keep the non-missing ones
    new_meta['country'] = new_meta['country_y'].mask(new_meta['country_y'].isna(), new_meta['country_x'])

    # Country: keep the non-missing ones
    new_meta['place'] = new_meta['place_y'].mask(new_meta['place_y'].isna(), new_meta['place_x'])

    # Clades: keep non-missing clades - subgenogroup
    new_meta['clade'] = new_meta['subgenogroup'].mask(new_meta['subgenogroup'].isna(), new_meta['clade_x'])

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
    new_meta['isolation_source'] = new_meta['isolation_source'].apply(standardize_isolation_source)
    new_meta['isolation'] = new_meta['isolation'].apply(standardize_isolation_source)

    # Combine the two columns, prioritizing the standardized values
    new_meta['combined_isolation_source'] = new_meta['isolation'].mask(new_meta['isolation'].isna(), new_meta['isolation_source'])

    # Replace the original isolation_source column
    new_meta['isolation_source'] = new_meta['combined_isolation_source']

    # Drop unnecessary columns
    new_meta = new_meta.drop(['isolation', 'combined_isolation_source'], axis=1)


    # Check if the regions file is supplied
    if con_reg_table:
        # Dictionary to store country-region mappings
        regions = {}
        
        # Read the regions file
        with open(con_reg_table) as f:
            regs = f.readlines()
        
        # Remove the first line (header)
        regs = regs[1:]
        
        # Populate the regions dictionary
        for x in regs:
            x = x.strip()
            pair = x.split("\t")
            if len(pair) > 1:
                regions[pair[0].strip().lower()] = pair[1].strip()
        
        # List to store the new region values
        newregion = []
        
        # Iterate over each country in the 'country' column of new_meta
        for coun in new_meta['country']:
            reg = "NA"  # Default value if no region is found
            if pd.notna(coun):
                coun = coun.strip().lower()
                coun_with_underscores = coun.replace(' ', '_')
                
                if coun not in regions:
                    if coun_with_underscores not in regions:
                        print(f"No region found for {coun}! Setting to NA")
                    else:
                        reg = regions[coun_with_underscores]
                else:
                    reg = regions[coun]
            
            reg = reg.replace('_', ' ')
            reg = reg.title()
            newregion.append(reg)
        
        # Update the 'region' column in the new_meta DataFrame with the new region values
        new_meta['region'] = newregion

    # Define a mapping for full terms to their abbreviations and standardized names
    short_versions = {
        'acute flaccid paralysis': 'AFP',
        'Hand-foot-and-mouth disease': 'HFMD',
        'central nervous system':'CNS',
        'Guillain-Barré syndrome': 'GBS',
        'Febrile Illness': 'Fever',
        'Febrile illness': 'Fever',
        'Cns symptoms': "CNS Symptoms",
        'Opsomyoclonus Syndrome': "OMS",
        'Myoclonic Jerk':'OMS',
        'Cns Involvement': "CNS Symptoms",
        'Cns Disorder': "CNS Symptoms",
        'Poliomyelitis-Like Disease': 'AFP',
        'Poliomyelitis-Like paralysis': 'AFP',
        'death': 'Fatality',
        'fatal': 'Fatality'
    }
    short_forms = set(short_versions.values())

    def clean_diagnosis(diagnosis, threshold=80):
        if pd.isna(diagnosis):
            return np.nan
        
        # Check if the diagnosis is already a short form
        if diagnosis in short_forms:
            return diagnosis
        
        # Remove punctuation and split multiple diagnoses
        clean_diag = diagnosis.replace(',', ';').replace('/', ';').replace('  ', ' ').strip()
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

    def extract_major_diagnosis(cleaned_diagnosis):
        if pd.isna(cleaned_diagnosis):
            return np.nan
        major_terms = {'AFP', 'HFMD', 'Encephalitis' ,'Fatality'}
        diagnoses = cleaned_diagnosis.split('; ')
        major_diagnoses = [diag for diag in diagnoses if diag in major_terms]
        return '; '.join(sorted(set(major_diagnoses)))

    # Apply the function to the 'Diagnosis' column
    new_meta['med_diagnosis_all'] = new_meta['diagnosis'].apply(lambda x: clean_diagnosis(x))

    # Apply the function to the 'Diagnosis' column
    new_meta['med_diagnosis_major'] = new_meta['med_diagnosis_all'].apply(lambda x: extract_major_diagnosis(x))
    
    # Add filter for age and add age ranges
    new_meta["has_age"] = ~new_meta["age_yrs"].isna()

    # parse gender
    new_meta['gender'] = new_meta['sex'].mask(new_meta['sex'].isna(),new_meta['gender'])
    new_meta['gender'] = new_meta['gender'].mask(new_meta['gender'].str.contains('female', case=False, na=False), 'F')
    new_meta['gender'] = new_meta['gender'].mask(new_meta['gender'].str.contains('male', case=False, na=False), 'M')

    #Define age bins and labels for years
    bins_years = [-np.inf, 1, 6, 11, 18, np.inf]
    labels_years = ['<1 y/o', '1-5 y/o', '6-10 y/o', '11-17 y/o', '18+ y/o']

    # Define age bins and labels for months
    bins_months = [-np.inf, 0.25, 0.5, 1]
    labels_months = ['0-3 m/o', '4-6 m/o', '7-12 m/o']

    # Create age_range column using pd.cut for years
    new_meta['age_range'] = pd.cut(new_meta['age_yrs'], bins=bins_years, labels=labels_years, right=False).astype(str)

    # Handle ages less than 1 year old separately using pd.cut for months
    mask_months = new_meta['age_yrs'] < 1
    new_meta.loc[mask_months, 'age_range'] = pd.cut(new_meta.loc[mask_months, 'age_yrs'], bins=bins_months, labels=labels_months, right=False).astype(str)

    # write new metadata file to output
    new_meta= new_meta.loc[:,['accession', 'genbank_accession_rev', 'strain', 'date', 'region', 'place',
        'country', 'host', 'gender', 'age_yrs','age_range', 'med_diagnosis_all','med_diagnosis_major',
        'isolation_source', 'length','date_submitted',
        'sra_accession', 'abbr_authors', 'reverse', 'authors', 'institution',
        'index','qc.overallScore', 'qc.overallStatus',
        'alignmentScore', 'alignmentStart', 'alignmentEnd', 'genome_coverage']]
    new_meta.to_csv(output_csv_meta, sep='\t', index=False)