import pandas as pd
import numpy as np
import re
import sys
import ipdb

#just adds very basic last-minute data (which cannot need reconstructing/ancestral) to the metadata
#before export for colouring/filtering etc
#ONLY SHOULD BE USED FOR TRAITS THAT DO NOT NEED RECONSTRUCTING!!

#Is purely to add metadata when something new comes in last-minute without re-running everything
#Format of the file:
#column 1: 'strain' or 'accession' to indicate how to indicate the record to be changed


#EX: to modify a strain called "SE02-21-01" so that the column 'symptom' has the value 'TM':
#strain	SE02-21-01	symptom	TM

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='add additional metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--add', help="tsv file with new data, format: 'column1 value1 column2 value2'. Column1 needs to be strain or accession")
    parser.add_argument('--input', help="input meta file")
    parser.add_argument('--rivm', help="input rivm subgenotype file")
    parser.add_argument('--id', help="id: strain or accession", choices=["strain","accession"],default="accession")
    parser.add_argument('--output', help="output meta file")
    args = parser.parse_args()

    new_data = pd.read_csv(args.add, sep='\t')
    meta = pd.read_csv(args.input, sep='\t', index_col=False)
    rivm_subtypes = pd.read_csv(args.rivm, index_col=False)
    id_field = args.id
    
    # ipdb.set_trace()
    

    if id_field not in new_data.columns:
            sys.exit(f"{id_field} was not in new dataframe.")
            
    # Remove duplicates based on id_field
    new_data = new_data.drop_duplicates(subset=id_field)
    meta = meta.drop_duplicates(subset=id_field)
    rivm_subtypes = rivm_subtypes.drop_duplicates(subset="name")
    
    # Select only the relevant columns from new_data
    new_data= new_data.loc[:,[id_field,"subgenogroup"]] # kept in case the RIVM subgenotypes are NA

    # Merge the dataframes on id_field
    new_meta = pd.merge(meta, new_data, 
                        on=id_field, 
                        how='left')

    # Add subgenotypes from RIVM
    mask = (rivm_subtypes["VP1 subgenogroup"]!="Could not assign")& (rivm_subtypes["type"]=="EV-A71")
    rivm_subtypes= rivm_subtypes.loc[mask, ["name","VP1 subgenogroup"]]

    # Change the colnames
    rivm_subtypes.rename(columns={"name":"accession", "VP1 subgenogroup":"RIVM_subgenogroup"}, inplace=True)

    # Merge the dataframes on id_field
    final_meta = pd.merge(new_meta, rivm_subtypes, on=id_field, how='left')

    # ipdb.set_trace()

    # save it to the new metadata file
    final_meta.to_csv(args.output, sep='\t', index=False)