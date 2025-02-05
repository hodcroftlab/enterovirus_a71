#!/bin/bash

# Input and output file paths from arguments
input_file="$1"
backup_file="$2"
output_file="$3"

temp_output_file="data/strain_names"

# Define the positions of accession and strain columns based on the header (1-based index)
accession_col=$(head -n 1 "$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++){if($i=="accession") print i}}')
strain_col=$(head -n 1 "$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++){if($i=="strain") print i}}')

# Check if the required columns are found
if [[ -z "$accession_col" || -z "$strain_col" ]]; then
    echo "Error: Could not find accession or strain column in the header."
    exit 1
fi

# Copy the backup file to the output file
if [[ -f "$backup_file" ]]; then
    cp "$backup_file" "$output_file"
    cp "$backup_file" "$temp_output_file"
else
    echo -e "accession\tstrain" > "$output_file"
fi

# Read existing accessions in the output file (from backup) into an array for faster lookup
existing_accessions=()
while IFS=$'\t' read -r accession strain; do
    existing_accessions["$accession"]=1
done < <(tail -n +2 "$output_file")  # Skip the header

# Cut out the accession and strain columns from input file and save to a temporary file
cut -f"$accession_col","$strain_col" "$input_file" | sort -u > "$temp_output_file"

# Check if the temp output file was created and is not empty
if [[ ! -s "$temp_output_file" ]]; then
    echo "Error: Temp output file $temp_output_file is empty or could not be created."
    exit 1
fi

skipped=0
new_accessions=0

# Process the temp file and append new entries to the output file
{
    # Skip the header line from the temp file
    read -r
    while IFS=$'\t' read -r accession strain; do
        # Skip lines where accession is empty
        if [[ -z "$accession" ]]; then
            continue
        fi

        # Check if accession already exists in the output file (using array lookup)
        if [[ ${existing_accessions["$accession"]} ]]; then
            ((skipped += 1))
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
                    ((new_accessions += 1))
                    echo "Updated strain for $accession: $strain"
                else
                    strain="$accession"
                fi
            else
                strain="$accession"
            fi
        fi

        # Append the updated entry to the output file
        echo -e "$accession\t$strain" >> "$output_file"
        existing_accessions["$accession"]=1
    done
} < "$temp_output_file"

# Sort and remove duplicates in the output file
sort -u "$output_file" -o "$output_file"

# Add the header at the end if not already present
if ! head -n 1 "$output_file" | grep -qP '^accession\s+strain'; then
    # Prepend the header to the output file
    echo -e "accession\tstrain" | cat - "$output_file" > temp && mv temp "$output_file"
fi

# Remove the last line from the output file
sed -i '$ d' "$output_file"

# Remove the temp file
rm "$temp_output_file"

# Report completion
echo -e "\nProcessing completed. Updated metadata saved to $output_file.\n$skipped accessions were skipped.\n$new_accessions new strains were added."