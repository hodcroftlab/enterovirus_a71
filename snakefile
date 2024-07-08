###############
# Snakemake execution templates:

# To run a default VP1 run(<600bp):
# snakemake  vp1/auspice/ev_a71_vp1.json --cores 1

# To run a default whole genome run ( <6400bp):
# snakemake whole_genome/auspice/ev_a71_whole_genome.json --cores 1

###############
# Define segments to analyze
segments = ['vp1', 'whole_genome']

# Expand augur JSON paths
augur_jsons = expand("auspice/ev_a71_{seg}.json", seg=segments)

rule all:
    input:
        augur_jsons

# Rule to handle configuration files
rule files:
    input:
        sequence_length =   "{seg}",
        dropped_strains =   "config/dropped_strains.txt",
        reference =         "{seg}/config/reference_sequence.gb",
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv"

files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

rule fetch:
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv"
    shell:
        """
            cd {input.dir} 
            snakemake --cores 5 all 
            cp {output.sequences} "../data/" 
            cp {output.metadata} "../data/" 
            cd ../
        """


#####################################################################################################
# BLAST
# blast fasta files for vp1 
rule blast:
    input: 
        blast_db_file = "data/references/reference_vp1_blast.fasta",
        seqs_to_blast = rules.fetch.output.sequences  # download from NCBI Virus
    output:
        blast_out = "vp1/temp/blast_out.csv"
    params:
        blast_db = "vp1/temp/entero_db_vp1"
    shell:
        """
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out, # output blast (vp1)
        input_seqs = rules.fetch.output.sequences  # download from NCBI Virus
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        matchLen = 300,
        range="{seg}"
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range}
        """

rule xlsx_to_csv:
    message:
        """
        Changing format from excel to csv. Clean up dates and other metadata.
        """
    input:
        dataframe =  "data/clade_assign_publications.xlsx"

    params:
        strain_id_field= "accession"
    output:
        dataframe = "data/clade_assign_publications.tsv"
    shell:
        """
        python scripts/xlsx2csv.py \
            --input {input.dataframe} \
            --id_field {params.strain_id_field}\
            --output {output.dataframe}
        """

rule add_metadata:
    message:
        """
        Cleaning dates in metadata
        """
    input:
        metadata =  "data/metadata.tsv",
        new_data =  rules.xlsx_to_csv.output.dataframe,
        regions = ancient(files.regions)

    params:
        strain_id_field= "accession"
    output:
        metadata = "{seg}/results/final_metadata.tsv"
    shell:
        """
        python scripts/add_metadata.py \
            --input {input.metadata} \
            --add {input.new_data} \
            --regions {input.regions}\
            --id {params.strain_id_field}\
            --output {output.metadata}
        """

# old way of cleaning up dates
rule curate: 
    message:
        """
        Cleaning dates in metadata
        """
    input:
        metadata = rules.add_metadata.output.metadata
    output:
        metadata = "{seg}/results/cleaned_metadata.tsv"
    shell:
        """
        python scripts/parse_date.py \
            --input {input.metadata} \
            --output {output.metadata}
        """

## use augur curate for formatting dates and others
# rule curate:
#     message:
#         """
#         Cleaning up metadata
#         """
#     input:
#         metadata =  rules.add_metadata.output.metadata
        
#     params:
#         strain_id_field= "accession",
#         date_column = ['date', 'date_submitted'],
#         format=['%Y', '%Y-%m', '%Y-%m-%d', '%Y-%m-%dT%H:%M:%SZ','%Y-XX-XX','%Y-%m-XX','XXXX-XX-XX']
#     output:
#         metadata = "{seg}/results/cleaned_metadata.tsv"
#     shell:
#         """
#         augur curate format-dates \
#             --metadata {input.metadata} \
#             --date-fields {params.date_column}\
#             --expected-date-formats {params.format}\
#             --id-column {params.strain_id_field}\
#             --output-metadata {output.metadata}
#         """


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.blast_sort.output.sequences
    output:
        sequence_index = "{seg}/results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.curate.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "{seg}/results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 2000,
        strain_id_field= "accession",
        min_date = 1965  # BrCr was collected in 1970
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --output {output.sequences}
        """

rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format
        """
    input:
        reference = files.reference

    output:
        reference = "{seg}/results/reference_sequence.fasta"
    shell:
        """
        python scripts/reference_genbank_to_fasta.py \
            --input {input.reference} \
            --output {output.reference}
        """

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "{seg}/results/aligned.fasta"

    params:
            nuc_mismatch_all = 10,
            nuc_seed_length = 30
    shell:
        """
        nextclade run \
        {input.sequences}  \
        --input-ref {input.reference}\
        --allowed-mismatches {params.nuc_mismatch_all} \
        --min-length {params.nuc_seed_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
        """

rule fix_align_codon:
    input:
        sequences = rules.align.output.alignment
    output:
        alignment = "{seg}/results/aligned_fixed.fasta"
    shell:
        """
        Rscript scripts/fixAlignmentGaps.R {input.sequences} {output.alignment}
        """

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.fix_align_codon.output.alignment
    output:
        tree = "{seg}/results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.fix_align_codon.output.alignment,
        metadata =  rules.curate.output.metadata,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 6, # was 3
        strain_id_field ="accession",
        clock_rate = 0.004,
        clock_std_dev = 0.0015
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        node_data = "{seg}/results/nt_muts.json",
        # reconstructed = "{seg}/results/nt_muts.json",
        # amino_acids = "{seg}/results/aa_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """
# --output-node-data {output.reconstructed} {output.amino_acids}\
# 
rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "{seg}/results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule clades: #check
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        clade_data = "{seg}/results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.curate.output.metadata
    output:
        node_data = "{seg}/results/traits.json",
    params:
        traits = "country",
        strain_id_field= "accession"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.traits} \
            --confidence
        """

rule clade_published:
    message: "Assigning clades from publications"
    input:
        metadata = rules.curate.output.metadata,
        new_data = "data/clade_assign_publications.tsv",
        rivm_data = "{seg}/data/{seg}_subgenotypes_rivm.csv"
    params:
        strain_id_field= "accession"
    output:
        meta = "{seg}/results/final_metadata_added_subgenotyp.tsv"
    shell:
        """
        python scripts/published_clades.py --input {input.metadata} --add {input.new_data} --rivm {input.rivm_data}\
        --id {params.strain_id_field} --output {output.meta}
        """

rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.clade_published.output.meta,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        # nt_muts = rules.ancestral.output.reconstructed,
        # aa_muts = rules.ancestral.output.amino_acids,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= "accession"
    output:
        auspice_json = "auspice/ev_a71_{seg}_accession.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule rename_json:
    input:
        auspice_json= rules.export.output.auspice_json,
        metadata = rules.clade_published.output.meta,
    output:
        auspice_json="auspice/ev_a71_{seg}.json"
    params:
        strain_id_field="accession",
        display_strain_field= "strain"
    shell:
        """
        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --metadata-id-columns {params.strain_id_field} \
                --input-auspice-json {input.auspice_json} \
                --display-strain-name {params.display_strain_field} \
                --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        # "results ",
        "auspice"
    shell:
        "rm -rfv {params}"


rule rename_whole_genome:
    message: "Removing directories: {params}"
    shell:
        """
        mv auspice/ev_a71_whole_genome.json auspice/ev_a71_whole-genome.json
        """
