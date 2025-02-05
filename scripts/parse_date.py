"""
date_parse.py standardizes dates contained in metadata.

Arguments:
    --input: the metadata.tsv file to clean
    --output: the cleaned metadata.tsv used for further analysis

date_parse.py is called within `snakefile`.

"""

import pandas as pd
import argparse

def standardize_date(oldDate):
    if pd.isna(oldDate):
        return "XXXX-XX-XX"
    elif oldDate.count("_") == 2:
        return oldDate.replace("_", "-")
    elif oldDate.count("_") == 1:
        return oldDate.replace("_", "-") + "-XX"
    elif oldDate.count("-") == 1:
        return oldDate + "-XX"
    elif len(oldDate) == 4:
        return oldDate + "-XX-XX"
    return oldDate

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='date_parse',
        description="Clean dates of metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar=' ', required=True, help="input metadata")
    parser.add_argument('-o', '--output', metavar=' ', required=True, help="output metadata")
    args = parser.parse_args()

    meta = pd.read_csv(args.input, keep_default_na=True, sep='\t', index_col=False)

    meta['date'] = meta['date'].apply(standardize_date)

    meta.to_csv(args.output, sep='\t', index=False)

