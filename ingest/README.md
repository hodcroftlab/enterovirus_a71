# Ingest

The `ingest` directory contains scripts and workflow files for data ingestion and preparation. Sequences 
and corresponding metdata are fetched from Genbank, then cleaned and curated. This process includes restructuring and formatting the metadata and sequence files into formats suitable for downstream processes.

## Directory Overview
- **[generate_from_genbank](bin/generate_from_genbank.py):** Script to download and parse GenBank files into required formats.
- **`Snakefile`:** Snakemake workflow file to run the ingest process.
- **`bin/`:** Directory containing scripts used in the ingest workflow.
- **`config`:** Configuration file for the ingest process.
- **`data/`:** Directory containing input data for the ingest process.
- **`source-data/`:** Directory for annotations and geo-location rules.
- **`vendored/`:** Directory for vendored scripts utilized in the ingest.
- **`workflow/`:** Directory containing the rules for the ingest workflow.



## First Steps

To set up and run the ingest workflow, follow these steps:

### 1. Prepare Reference Files
[Nextclade](https://clades.nextstrain.org/) is used as part of the `ingest` workflow to align sequences to the reference, and assign the sequences into clades. Nextclade requires reference files to be in a specific format, such as `reference.fasta` and `annotation.gff3`. The following steps can be used to generate these specific file formats from the reference sequence which is in a Genbank file format.

1. **Verify `config` Settings:**  
   Open the `ingest/config/config.yaml` file and confirm the `ncbi_taxon_id` is correct.

2. **Run `generate_from_genbank.py` Script:**  
   Execute the script (located in ingest/bin/) to generate required reference files:
   ```bash
   python3 ingest/generate_from_genbank.py --reference "U22521.1" --output-dir ingest/data/references/
   ```

During execution, you may be provided to provide CDS annotations. You can specify one of the following fields to use as the CDS names:

- `[0]`
- `[product]` or `[leave empty for manual choice]` to select proteins.
- `[2]`.

Please note that for some pathogens, no common fields may be available in the set. In such cases, you will need to leave the input blank and manually assign CDS names.  
The generated files will be saved in the `data/references` subdirectory and used by the `ingest` Snakefile.

3. **Update Attributes**  
   Ensure that attributes in `data/references/pathogen.json` are up-to-date. Please consult the [Nextclade pathogen configuration documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html)  for more optional attributes to add to the file.  



## 2. Run the Ingest Workflow

Run the ingest workflow using the Snakefile. Depending on your system, you may need to make the scripts executable first:

```bash
chmod +x ./vendored/*; chmod +x ./bin/*
```
The Snakefile may be run by executing the following command in your terminal from within the ingest directory:

```bash
snakemake --cores 1 all
```


## Updating Vendored Scripts

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage vendored scripts in `ingest/vendored`.

### Steps to Update Vendored Scripts

1. Install `git subrepo` by following the [installation guide](https://github.com/ingydotnet/git-subrepo#installation).
2. Pull the latest changes from the central ingest repository by following the instructions in [`ingest/vendored/README.md`](vendored/README.md#vendoring).

