# Enterovirus A71 Nextstrain Analysis

This repository provides a comprehensive Nextstrain analysis of Enterovirus A71. You can choose to perform either a **VP1 run (>=600 base pairs)** or a **whole genome run (>=6400 base pairs)**.

For those unfamiliar with Nextstrain or needing installation guidance, please refer to the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/).

### Enhancing the Analysis
This analysis would benefit from additional metadata, such as patient age, spatial data, and clinical outcomes. If you have relevant data and are willing to share, please [contact us](mailto:nadia.neuner-jehle@swisstph.ch).

The data for this analysis is available from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). Instructions for downloading sequences are provided under [Sequences](#sequences).

## Repository Organization
This repository includes the following directories and files:

- `scripts`: Custom Python scripts called by the `snakefile`.
- `snakefile`: The entire computational pipeline, managed using Snakemake. Snakemake documentation can be found [here](https://snakemake.readthedocs.io/en/stable/).
- `ingest`: Contains Python scripts and the `snakefile` for automatic downloading of EV-A71 sequences and metadata.
- `vp1`: Sequences and configuration files for the **VP1 run**.
- `whole_genome`: Sequences and configuration files for the **whole genome run**.

### Configuration Files
The `config`, `vp1/config`, and `whole_genome/config` directories contain necessary configuration files:
- `colors.tsv`: Color scheme
- `geo_regions.tsv`: Geographical locations
- `lat_longs.tsv`: Latitude data
- `dropped_strains.txt`: Dropped strains
- `clades_genome.tsv`: Virus clade assignments
- `reference_sequence.gb`: Reference sequence
- `auspice_config.json`: Auspice configuration file

The reference sequence used is [BrCr, accession number U22521](https://www.genome.jp/dbget-bin/www_bget?genbank-vrl:U22521), sampled in 1970.

## Quickstart

### Setup

#### Nextstrain Environment
Install the Nextstrain environment by following [these instructions](https://docs.nextstrain.org/en/latest/guides/install/local-installation.html).

### Running a Build

Activate the Nextstrain environment:
```bash
micromamba activate nextstrain
```

To perform a build, run:
```bash
snakemake --cores 9 all
```

For specific builds:
- VP1 build:
```bash
snakemake auspice/ev_a71_vp1.json --cores 9
```
- Whole genome build:
```bash
snakemake auspice/ev_a71_whole-genome.json --cores 9
```

For tanglegrams, we can run the build on sub-alignments of the whole genome alignment. 
You can either run it for the specific genes or for the proteins P1, P2, P3.
- gene build:
```bash
snakemake all_genes --cores 9
```
- Whole genome build:
```bash
snakemake all_proteins --cores 9
```

> [!NOTE]
> Version of <ins> augur</ins>: `augur 27.0.0`\
> Version of <ins> auspice</ins>: `auspice 2.59.1`

## Ingest
For more information on how to run the `ingest`, please refer to the [README](ingest/README.md) in the `ingest` folder.

### Visualizing the Build
To visualize the build, use Auspice:
```bash
auspice view --datasetDir auspice
```
To run two visualizations simultaneously, you may need to set the port:
```bash
export PORT=4001
```

### Sequences
Sequences can be downloaded manually or automatically.

1. **Manual Download**: Visit [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), search for `EV-A71` or Taxid `39054`, and download the sequences.
2. **Automated Download**: The `ingest` functionality, included in the main `snakefile`, handles automatic downloading.

The ingest pipeline is based on the Nextstrain [RSV ingest workflow](https://github.com/nextstrain/rsv.git). Running the **ingest** pipeline produces `data/metadata.tsv` and `data/sequences.fasta`.

## Feedback
For questions or comments, contact me via GitHub or [nadia.neuner-jehle@swisstph.ch](mailto:nadia.neuner-jehle@swisstph.ch)

## To Do:
- [X] Overwrite NCBI virus metadata with "corrected" collection dates
- [X] Replace `parse_date` with `augur curate`
- [X] Provide a way to create and use "local" accession numbers for sequences not on Genbank yet.
- [ ] Update symptom list
- [ ] Get [update_strain.sh](scripts/update_strain.sh) to work


## Acknowledgments
- [Nextstrain](https://nextstrain.org/)
- [Auspice](https://auspice.us/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Biopython](https://biopython.org/)
- [Genbank](https://www.ncbi.nlm.nih.gov/genbank/)
- [NCBI](https://www.ncbi.nlm.nih.gov/)