# Ingest

The `ingest` directory contains scripts and workflow files for data ingestion and preparation.

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
You will need specific reference files, such as `reference.fasta` and `annotation.gff3`.

1. **Verify `config` Settings:**  
   Open the `config` file and confirm the `taxid` is correct.

2. **Run `generate_from_genbank.py` Script:**  
   Execute the script to generate required reference files:
   ```bash
   python3 ingest/generate_from_genbank.py --reference "U22521.1" --output-dir ingest/data/references/
   ```

During execution, you may need to specify the following:

- `[0]`
- `[product]` or `[leave empty for manual choice]` to select proteins.
- `[2]`.

The generated files will be saved in the `data/references` subdirectory and used by the `ingest` Snakefile.

3. **Update Attributes**  
   Ensure that attributes in `data/references/pathogen.json` are up-to-date.



## 2. Run the Ingest Workflow

Run the ingest workflow using the Snakefile. Depending on your system, you may need to make the scripts executable first:

```bash
chmod +x ./vendored/*; chmod +x ./bin/*
```

## Updating Vendored Scripts

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage vendored scripts in `ingest/vendored`.

### Steps to Update Vendored Scripts

1. Install `git subrepo` by following the [installation guide](https://github.com/ingydotnet/git-subrepo#installation).
2. Pull the latest changes from the central ingest repository by following the instructions in [`ingest/vendored/README.md`](vendored/README.md#vendoring).

