#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from Bio import SeqIO
# import ipdb


def normalize_strain(name):
    """Normalize a strain name so that trivial formatting differences
    (e.g. '_' vs '|', extra whitespace, mixed case) don't stop us from
    recognizing two records as the same strain."""
    if pd.isna(name):
        return name
    name = str(name).strip().upper()
    # Collapse common separator variants ('_', '|', '/', '-', '.', whitespace)
    # down to a single canonical separator so 'CF210042_FRA06' and
    # 'CF210042|FRA06' are treated as the same strain.
    name = re.sub(r"[_|/.\s]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name


def extract_year(date_str):
    """Pull a 4-digit year out of a date string, if present."""
    if pd.isna(date_str):
        return None
    match = re.search(r"(1[89]\d{2}|20\d{2})", str(date_str))
    return match.group(1) if match else None


def normalize_country(country):
    if pd.isna(country) or not str(country).strip():
        return None
    return str(country).strip().upper()


def percent_identity(seq_a, seq_b):
    """Simple ungapped percent identity over the aligned (shorter) length,
    normalized by the longer sequence's length so that a shorter sequence
    fully contained in a longer one still scores lower if the longer one
    has extra unmatched content."""
    length = min(len(seq_a), len(seq_b))
    if length == 0:
        return 0.0
    matches = sum(1 for a, b in zip(seq_a[:length], seq_b[:length]) if a == b)
    return matches / max(len(seq_a), len(seq_b))


def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate sequences with identical (or near-identical) strain names AND identical sequences"
    )
    parser.add_argument("--in-sequences", required=True, help="Input FASTA")
    parser.add_argument("--metadata", required=True, help="Input metadata TSV")
    parser.add_argument("--out-sequences", required=True, help="Output FASTA")
    parser.add_argument("--id-field", default="accession", help="Metadata column matching FASTA record IDs (default: accession)")
    parser.add_argument("--strain-field", default="strain", help="Metadata column containing strain names for deduplication (default: strain)")
    parser.add_argument("--country-field", default="country", help="Metadata column containing country (default: country)")
    parser.add_argument("--date-field", default="date", help="Metadata column containing collection date (default: date)")
    parser.add_argument("--threshold", type=float, default=0.995, help="Percent identity threshold to consider sequences identical when country/year cannot be compared (default: 0.995)")
    parser.add_argument("--relaxed-threshold", type=float, default=0.95, help="Lower percent identity threshold to use when country AND year both match between the two records (default: 0.95)")
    parser.add_argument("--log-file", default=None, help="Optional TSV path to log removed duplicate accessions and why they were removed")

    args = parser.parse_args()
    # ------------------
    # Read metadata
    # ------------------
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    if args.id_field not in meta.columns:
        raise ValueError(
            f"ID field '{args.id_field}' not found in metadata columns"
        )

    # ------------------
    # Read FASTA
    # ------------------
    records = list(SeqIO.parse(args.in_sequences, "fasta"))

    fasta_df = pd.DataFrame(
        {
            args.id_field: [rec.id for rec in records],
            "sequence": [str(rec.seq) for rec in records],
        }
    )

    # Metadata accession, strain, country, date, length
    cols = [args.id_field, args.strain_field, "NCBI_length_genome"]
    optional_cols = {"country": args.country_field, "date": args.date_field}
    for out_name, col in optional_cols.items():
        if col in meta.columns:
            cols.append(col)

    meta_df = meta[cols].copy()
    rename_map = {args.id_field: "accession", args.strain_field: "strain", "NCBI_length_genome": "length"}
    for out_name, col in optional_cols.items():
        if col in meta.columns:
            rename_map[col] = out_name
    meta_df = meta_df.rename(columns=rename_map)

    meta_df = meta_df.loc[meta_df["accession"].isin(fasta_df[args.id_field])]
    meta_df["length"] = pd.to_numeric(meta_df["length"], errors="coerce")

    # Normalize strain names so that separator differences (e.g. '_' vs '|')
    # don't hide duplicates from each other.
    meta_df["strain_norm"] = meta_df["strain"].apply(normalize_strain)

    if "country" in meta_df.columns:
        meta_df["country_norm"] = meta_df["country"].apply(normalize_country)
    else:
        meta_df["country_norm"] = None

    if "date" in meta_df.columns:
        meta_df["year"] = meta_df["date"].apply(extract_year)
    else:
        meta_df["year"] = None

    # Sort by length descending so the longest sequence in each strain
    # group is always the first (reference/kept) record.
    meta_df = meta_df.sort_values("length", ascending=False)

    # Groups of accessions that share a normalized strain name, including
    # the longest ("reference") record -- not just the extra duplicates.
    strain_counts = meta_df["strain_norm"].value_counts()
    dup_strain_names = strain_counts[strain_counts > 1].index

    acc_arrays = (
        meta_df[meta_df["strain_norm"].isin(dup_strain_names)]
        .groupby("strain_norm", sort=False)["accession"]
        .apply(list)
        .reset_index(name="accession_array")
    )

    meta_by_acc = meta_df.set_index("accession")

    sequences_to_remove = set()
    removal_log = []
    # Percent identity of sequences with the same (normalized) strain name
    for strain, acc in acc_arrays.itertuples(index=False):
        seqs_by_acc = fasta_df.set_index(args.id_field)["sequence"]
        ref_acc = acc[0]
        ref_seq = seqs_by_acc.get(ref_acc)
        ref_country = meta_by_acc.loc[ref_acc, "country_norm"]
        ref_year = meta_by_acc.loc[ref_acc, "year"]

        for other_acc in acc[1:]:
            seq = seqs_by_acc.get(other_acc)
            if ref_seq is None or seq is None:
                continue

            pid = percent_identity(ref_seq, seq)

            # If both country and year are known and agree between the two
            # records, we're more confident they're the same isolate, so
            # allow a lower percent-identity threshold to still call it a
            # duplicate (accounting for sequencing/assembly noise).
            other_country = meta_by_acc.loc[other_acc, "country_norm"]
            other_year = meta_by_acc.loc[other_acc, "year"]
            same_country_year = (
                ref_country is not None
                and other_country is not None
                and ref_country == other_country
                and ref_year is not None
                and other_year is not None
                and ref_year == other_year
            )
            effective_threshold = args.relaxed_threshold if same_country_year else args.threshold

            if pid >= effective_threshold:
                sequences_to_remove.add(other_acc)  # keep the longest (ref_acc)
                removal_log.append(
                    {
                        "strain": strain,
                        "kept_accession": ref_acc,
                        "kept_length": len(ref_seq),
                        "removed_accession": other_acc,
                        "removed_length": len(seq),
                        "percent_identity": round(pid * 100, 2),
                        "threshold_used": effective_threshold,
                        "same_country_year": same_country_year,
                    }
                )

    # ------------------
    # Deduplicate
    # ------------------
    dedup_fasta_df = fasta_df[
        ~fasta_df[args.id_field].isin(sequences_to_remove)
    ]

    keep_ids = set(dedup_fasta_df[args.id_field])

    with open(args.out_sequences, "w") as out_fasta:
        SeqIO.write(
            (rec for rec in records if rec.id in keep_ids),
            out_fasta,
            "fasta",
        )

    # ------------------
    # Log removed duplicates
    # ------------------
    if args.log_file:
        log_df = pd.DataFrame(
            removal_log,
            columns=[
                "strain",
                "kept_accession",
                "kept_length",
                "removed_accession",
                "removed_length",
                "percent_identity",
                "threshold_used",
                "same_country_year",
            ],
        )
        log_df.to_csv(args.log_file, sep="\t", index=False)

    # ------------------
    # Summary
    # ------------------
    print(f"Input sequences: {len(records)}")
    print(f"Deduplicated sequences: {len(dedup_fasta_df)}")
    print(f"Removed duplicates: {len(records) - len(dedup_fasta_df)}")
    if args.log_file:
        print(f"Removed duplicate log written to: {args.log_file}")


if __name__ == "__main__":
    main()
