#!/bin/bash

# Input and output file paths from arguments
input_file="$1"
temp_output_file="data/strain_names"
output_file="$2"

# Define the positions of accession and strain columns based on the header (1-based index)
accession_col=$(head -n 1 "$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++){if($i=="accession") print i}}')
strain_col=$(head -n 1 "$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++){if($i=="strain") print i}}')

# Check if the required columns are found
if [[ -z "$accession_col" || -z "$strain_col" ]]; then
    echo "Error: Could not find accession or strain column in the header."
    exit 1
fi

# Cut out the accession and strain columns and save to a temporary file
cut -f"$accession_col","$strain_col" "$input_file" | sort -u > "$temp_output_file"

# Check if the temp output file was created and is not empty
if [[ ! -s "$temp_output_file" ]]; then
    echo "Error: Temp output file $temp_output_file is empty or could not be created."
    exit 1
fi

# Create or update the output file
if [[ ! -f "$output_file" ]]; then
    echo -e "accession\tstrain" > "$output_file"
fi

cat "$temp_output_file" >> "$output_file"
sort -u "$output_file" -o "$output_file"

# Load existing accessions from the output file
declare -A existing_strains
while IFS=$'\t' read -r accession strain; do
    existing_strains["$accession"]="$strain"
done < <(tail -n +2 "$output_file")

# Process the temp file and append new entries to the output file
{
    # Skip the header line from the temp file
    read -r
    while IFS=$'\t' read -r accession strain; do
        # Skip lines where accession is empty
        if [[ -z "$accession" ]]; then
            continue
        fi

        # Update the strain if needed
        if [ "$accession" == "$strain" ]; then
            # Fetch the GenBank record
            gb_entry=$(efetch -db nucleotide -id "$accession" -format gb 2>/dev/null)

            if [ $? -eq 0 ]; then
                # Extract the strain from the FEATURES section
                strain_from_gb=$(echo "$gb_entry" | grep -oP '/strain="\K[^"]+')

                if [ -n "$strain_from_gb" ]; then
                    strain=$strain_from_gb
                    echo "Updated strain for $accession: $strain"
                else
                    strain="$accession"
                fi
            else
                strain="$accession"
            fi

            # Update the output file
            echo -e "$accession\t$strain" >> "$output_file"
            sort -u "$output_file" -o "$output_file"
        fi
    done
} < "$temp_output_file"

# Sort and remove duplicates in the output file
sort -u "$output_file" -o "$output_file"

# Add the header at the end if not already present
if ! head -n 1 "$output_file" | grep -qP '^accession\s+strain'; then
    # Prepend the header to the output file
    echo -e "accession\tstrain" | cat - "$output_file" > temp && mv temp "$output_file"
fi

# Remove the temp file
rm "$temp_output_file"

# Report completion
echo "Processing completed. Updated metadata saved to $output_file"
