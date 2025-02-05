import pandas as pd
import numpy as np
import re
import ipdb
import sys

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='format additional metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  metavar=' ', help="input xslx")
    parser.add_argument('--id_field', help="strain or accession", choices=["strain","accession"],default="accession")
    parser.add_argument('-o', '--output', metavar=' ', help="output csv")
    args = parser.parse_args()

    input_xlsx = args.input
    id = args.id_field
    output_csv = args.output

    
    data_csv =pd.read_excel(io=input_xlsx)
    data_csv["Collection yr"]=data_csv["Collection yr"].astype('Int64')
    
    # ipdb.set_trace()
    # Assuming 'data_csv' is your DataFrame
    # data_csv = pd.DataFrame(...)

    # List of current column names from the DataFrame
    column_names = data_csv.columns.tolist()

    # Cleaned column names
    cleaned_column_names = []

    for col in column_names:
        # Convert to lowercase
        col = col.lower()
        # Remove parentheses
        col = col.replace('(', '').replace(')', '')
        # Replace spaces with underscores
        col = col.replace(' ', '_')
        # Remove any special characters except for underscores
        col = re.sub(r'\W+', '', col)
        cleaned_column_names.append(col)

    # Replace the data_csvsting column names with the new ones
    data_csv.columns = cleaned_column_names

    # Accession_no is special, replace that seperately
    data_csv.rename(columns={"accession_no":"accession"},inplace=True)

    # Process the collection_date column
    data_csv.collection_date=data_csv.collection_date.map(
        lambda date_str: f"{date_str.split('-')[2]}-XX-XX" if isinstance(date_str, str) and 'XX' in date_str else (
            f"{date_str.split('-')[2]}/{date_str.split('-')[1]}/{date_str.split('-')[0]}" if isinstance(date_str, str) else date_str
        )
    )

    
    # Identify duplicated rows based on the "accession" column
    duplicated_accessions = data_csv[data_csv.duplicated(id, keep=False)]

    # Function to count non-null values
    def count_non_null(row):
        return row.count()

    # Group by 'accession' and keep the row with the most non-null values
    best_rows = duplicated_accessions.groupby(id).apply(lambda group: group.loc[group.apply(count_non_null, axis=1).idxmax()])
    
    # Remove the index added by the groupby operation
    best_rows.reset_index(drop=True, inplace=True)

    # Drop duplicated rows from the original DataFrame (data_csv) and append the best rows
    result_data_csv = pd.concat([data_csv.drop_duplicates(id, keep=False), best_rows]).sort_index()
    # result_data_csv.head()

    result_data_csv.to_csv(output_csv, sep="\t", index=None)