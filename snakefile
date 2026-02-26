###############
# Snakemake execution templates:

# To run a default VP1 run:
# snakemake  auspice/enterovirus_A71_vp1.json --cores 9

# To run a default whole genome run (>6400bp):
# snakemake auspice/enterovirus_A71_whole-genome.json --cores 9

# To run specific genes for tanglegrams:
# snakemake --cores 9 all_genes

# To run specific proteins for tanglegrams:
# snakemake --cores 9 all_proteins

###############
TAXID = "39054"  # NCBI Taxonomy ID for Enterovirus A71
REFERENCE = "U22521"

if not config:
    configfile: "config/config.yaml"

import os
from datetime import date

# Load environment variables
# Try to load .env, but don't fail if it doesn't exist (for Actions)
try:
    from dotenv import load_dotenv
    load_dotenv(".env")
except:
    pass

REMOTE_GROUP = os.getenv("REMOTE_GROUP")
UPLOAD_DATE = date.today().isoformat()

###############
wildcard_constraints:
    seg="vp1|whole_genome|P1",
    gene="|-5utr|-vp4|-vp2|-vp3|-vp1|-2A|-2B|-2C|-3A|-3B|-3C|-3D|-3utr",
    protein="|-P1|-P2|-P3"
   
# Define segments to analyze
segments = ['vp1', 'whole-genome', "P1"]
GENES=["-5utr","-vp4", "-vp2", "-vp3", "-vp1", "-2A", "-2B", "-2C", "-3A", "-3B", "-3C", "-3D","-3utr"]
PROT = ["-P1", "-P2", "-P3"]
CODING_GENES = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
PROTEIN1 = ["VP4", "VP2", "VP3", "VP1"]

# parameters
INCL_ENPEN = True # include or exclude ENPEN data
FETCH_SEQUENCES = True
INCLUDE_REF = False

# Rule to handle configuration files
rule files:
    input:
        dropped_strains =   "config/exclude.txt",
        incl_strains =      "config/include.txt",   
        reference =         "config/reference_sequence.gb",
        gff_reference =     "{seg}/config/annotation.gff3",
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv",
        
        meta_public=        "data/meta_public.tsv",
        meta_collab =       "data/meta_ENPEN.tsv",
        meta_genbank =      "data/genbank_metadata.tsv",
        last_updated_file = "data/date_last_updated.txt",
        local_accn_file =   "data/local_accn.txt",
        SEQUENCES =         "data/sequences.fasta",
        METADATA =          "data/metadata.tsv",
        RIVM_CLADES =       "data/subgenotypes_rivm.csv",
        VP1_CLADES =        "data/clades_vp1.tsv",



files = rules.files.input

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/enterovirus_A71_{segs}.json", segs=segments),
        epitopes = "vp1/results/epitopes.json",
        meta = files.METADATA,
        seq = files.SEQUENCES


rule all_genes:
    input:
        seg_jsn = expand("auspice/enterovirus_A71_{segs}.json", segs=segments),
        augur_jsons = expand("auspice/enterovirus_A71_gene_{genes}.json", genes=GENES),
        epitopes = "vp1/results/epitopes.json",
        meta = files.METADATA,
        seq = files.SEQUENCES
        

rule all_proteins:
    input:
        augur_jsons = expand("auspice/enterovirus_A71_protein_{proteins}.json", proteins=PROT),
        seg_jsn = expand("auspice/enterovirus_A71_{segs}.json", segs=segments),
        epitopes = "vp1/results/epitopes.json",
        meta = files.METADATA,
        seq = files.SEQUENCES

rule next_update:
    input:
        "vp1/results/epitopes.json",
        "auspice/enterovirus_A71_vp1.json", 
        "auspice/enterovirus_A71_whole-genome.json"
        # "auspice/enterovirus_A71_gene_-vp1.json", 
        # "auspice/enterovirus_A71_gene_-3D.json"

rule ENPEN:
    input:
        "auspice/enterovirus_A71_P1.json", 


##############################
# Download from NBCI Virus with ingest snakefile
###############################
if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=files.SEQUENCES,
            metadata=files.METADATA
        threads: workflow.cores
        shell:
            """
            cd {input.dir} 
            snakemake --cores {threads} all
            cd ../
            """


# This rule is very slow. Only give accessions as input where you are certain that they have GenBank metadata.
rule fetch_metadata:
    message:
        """
        Retrieving GenBank metadata for the specified accessions. See {log} for details.
        """
    input:
        accessions="data/metadata/genbank_afm_failed.txt",
        config="config/config.yaml", # include symptom list and isolation source mapping
        lat_longs=files.lat_longs
    output:
        metadata="data/metadata/genbank_afm_cases.tsv",
    params:
        virus="Enterovirus A71",
        genbank_metadata=files.meta_genbank,
        cols = ["strain", "accession", "country", "place", "region", "subgenogroup", "lineage", "date", "collection_yr", "gender", "age_yrs", "age_mo", "diagnosis", "isolation", "origin", "doi"],
    log:
        "logs/fetch_metadata.log"
    shell:
        """
        python scripts/fetch_genbank_metadata.py \
            --virus "{params.virus}" \
            --accession_file {input.accessions} \
            --output {output.metadata} \
            --genbank {params.genbank_metadata} \
            --config {input.config} \
            --latlongs {input.lat_longs} \
            --columns {params.cols} \
            2> {log}
        """

##############################
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
###############################

rule curate:
    message:
        """
        Cleaning up metadata with augur curate
        """
    input:
        metadata=files.meta_public,  # Path to input metadata file
        meta_collab = files.meta_collab,  # Data shared with us by collaborators
        meta_genbank = files.meta_genbank
    params:
        strain_id_field=config["id_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
    output:
        merge = temp("data/merge_meta.tsv"),  # Final output file for publications metadata
        meta="data/curated/all_meta.tsv"  # Final merged output file
    shell:
        """        
        # Merge curated metadata
        augur merge --metadata metadata={input.metadata} meta_collab={input.meta_collab} meta_genbank={input.meta_genbank} \
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.merge}
        
        # Normalize strings and format dates for metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {output.merge} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.meta}
        echo "Curated metadata saved to {output.meta}"
        """


##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################

rule update_sequences:
    input:
        sequences = files.SEQUENCES,
        metadata=files.METADATA,
        extra_metadata = rules.curate.output.meta
    output:
        sequences = "data/all_sequences.fasta"
    params:
        file_ending = "data/*.fas*",
        temp = "data/temp_sequences.fasta",
        date_last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    shell:
        """
        touch {params.temp} && rm {params.temp}
        cat {params.file_ending} > {params.temp}
        python scripts/update_sequences.py --in_seq {params.temp} --out_seq {output.sequences} --dates {params.date_last_updated} \
        --local_accession {params.local_accn} --meta {input.metadata} --add {input.extra_metadata} \
        --ingest_seqs {input.sequences}
        rm {params.temp}
        awk '/^>/{{if (seen[$1]++ == 0) print; next}} !/^>/{{print}}' {output.sequences} > {params.temp} && mv {params.temp} {output.sequences}
        """



##############################
# BLAST
# blast fasta files for your specific proteins
# cut out your protein from fasta sequences
###############################

rule extract:
    input: 
        genbank_file = files.reference
    output: 
        extracted_fasta = "{seg}/config/reference.fasta",    
        extracted_genbank = "{seg}/config/reference.gbk",
    params:
        product_name = "{seg}",
        taxid = TAXID,
        annotation = lambda wildcards: f'--output_gff {wildcards.seg}/config/annotation.gff3' if wildcards.seg != "whole_genome" else ""

    shell:
        """
        python scripts/extract_gene_from_whole_genome.py \
        --genbank_file {input.genbank_file} \
        --output_fasta {output.extracted_fasta} \
        --product_name {params.product_name} \
        --output_genbank {output.extracted_genbank} \
        --taxid {params.taxid} \
        {params.annotation}
        """

rule blast:
    input: 
        blast_db_file = rules.extract.output.extracted_fasta,  
        seqs_to_blast = rules.fetch.output.sequences
    output:
        blast_out = "temp/{seg}/blast_out.csv"
    params:
        blast_db =  "temp/{seg}/blast_database"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} \
        -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out, # output blast (for your protein)
        input_seqs = rules.fetch.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta",
    params:
        range = "{seg}",  # Determines which protein (or whole genome) is processed
        min_length = lambda wildcards: {"vp1": 600, "whole_genome": 6400, "P1": 2000}[wildcards.seg],  # Min length
        max_length = lambda wildcards: {"vp1": 900, "whole_genome": 8000, "P1": 2600}[wildcards.seg]  # Max length
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range} \
            --min_length {params.min_length} \
            --max_length {params.max_length}
        """


##############################
# Merge all metadata files (NCBI download and own files) and clean them up
# potentially use augur merge: but not the same output can be achieved with augur
###############################

rule add_metadata:
    message:
        """
        Cleaning data in metadata
        """
    input:
        metadata=files.METADATA,
        new_data=rules.curate.output.meta,
        regions=ancient(files.regions),
        C1like_accn = "data/list_C1like_accn.txt",
    params:
        strain_id_field=config["id_field"],
        last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    output:
        metadata="data/all_metadata.tsv"
    shell:
        """
        python scripts/add_metadata.py \
            --input {input.metadata} \
            --add {input.new_data} \
            --local {params.local_accn} \
            --C1like {input.C1like_accn} \
            --update {params.last_updated}\
            --regions {input.regions} \
            --id {params.strain_id_field} \
            --output {output.metadata}
        """

## Deduplicate sequences that have identical strain names and sequences
rule deduplicate:
    message:
        """
        Segment: {wildcards.seg}
        Deduplicating sequences with identical strain names and sequences
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        metadata = rules.add_metadata.output.metadata
    params:
        id_field = config["id_field"],
        threshold = 0.995 # percent identity threshold to consider sequences as duplicates
    output:
        sequences = "{seg}/results/deduplicated_sequences.fasta",
    shell:
        """
        python scripts/deduplicate.py \
            --in-sequences {input.sequences} \
            --metadata {input.metadata} \
            --id-field {params.id_field} \
            --threshold {params.threshold} \
            --out-sequences {output.sequences} 
        """


##############################
# Create an index of sequence composition for filtering & filter
###############################
rule index_sequences:
    message:
        """
        Segment: {wildcards.seg}
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.deduplicate.output.sequences
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
        Segment: {wildcards.seg}

        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.deduplicate.output.sequences, ## x had no sequence data -> dropped because they don't meet the min sequence length in blast_sort
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.add_metadata.output.metadata,
        exclude = files.dropped_strains,
        include = files.incl_strains,
    output:
        sequences = "{seg}/results/filtered.fasta",
        reason ="{seg}/results/reasons.tsv",
        include = "{seg}/results/include.txt"
    log:
        "logs/filter.{seg}.log"
    params:
        group_by = "country year", # dropped because of ambiguous year information
        sequences_per_group = 500, # set lower if you want to have a max sequences per group
        strain_id_field = config["id_field"],
        min_date = 1975,  # BrCr was collected in 1970
        exclude_enpen = "--exclude-where ENPEN=True" if not INCL_ENPEN else "", # INCL_ENPEN was defined on line 32
        ref = REFERENCE if INCLUDE_REF else "" # add reference sequence to the include list if you want to keep it in the tree
    shell:
        """
        cat {input.include} > {output.include}
        echo "{params.ref}" >> {output.include}
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            {params.exclude_enpen} \
            --include {output.include} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --output-sequences {output.sequences}\
            --output-log {output.reason} \
            >> {log} 2>&1

        echo "Filtered sequences saved to {output.sequences}"
        """

##############################
# Reference sequence &
# Alignment
###############################
rule align: 
    message:
        """
        Segment: {wildcards.seg}
        Aligning sequences to {input.reference} using Nextclade.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = rules.filter.output.sequences,
        reference = rules.extract.output.extracted_fasta,
    output:
        alignment = "{seg}/results/aligned.fasta",
        tsv = "{seg}/results/nextclade.tsv",    
    benchmark:
        "benchmark/align.{seg}.log"
    params:
        penalty_gap_extend = config["align"]["penalty_gap_extend"],
        penalty_gap_open = config["align"]["penalty_gap_open"],
        penalty_gap_open_in_frame = config["align"]["penalty_gap_open_in_frame"],
        penalty_gap_open_out_of_frame = config["align"]["penalty_gap_open_out_of_frame"],
        kmer_length = config["align"]["kmer_length"],
        kmer_distance = config["align"]["kmer_distance"],
        min_match_length = config["align"]["min_match_length"],
        allowed_mismatches = config["align"]["allowed_mismatches"],
        min_length = config["align"]["min_length"]
        ## min_length
    threads: workflow.cores
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations "{wildcards.seg}/results/translations/cds_{{cds}}.translation.fasta" \
        --output-fasta {output.alignment}
        """

#  one-by-one genes
rule sub_alignments:
    message:
        """
        Segment: {wildcards.seg}
        Creating sub-alignments for {wildcards.gene}{wildcards.protein}
        """
    input:
        alignment=rules.align.output.alignment,
        reference=files.reference
    output:
        alignment = "{seg}/results/aligned{gene}{protein}.fasta"
    benchmark:
        "benchmark/sub_alignments.{seg}{gene}{protein}.log"
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq

        if wildcards.protein:
            real_gene = wildcards.protein.replace("-", "", 1)
            boundaries = {
                'P1': (744, 3329),
                'P2': (3330, 5063),
                'P3': (5064, 7322)
            }
            b = boundaries[real_gene]
        else:
            real_gene = wildcards.gene.replace("-", "", 1)

            # Extract boundaries from the reference GenBank file
            gene_boundaries = {}
            with open(input.reference) as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS" and 'Name' in feature.qualifiers:
                            product = feature.qualifiers['Name'][0].upper()
                            if product == real_gene.upper():
                                # Corrected: Use .start and .end directly
                                gene_boundaries[product] = (feature.location.start, feature.location.end)

            if real_gene.upper() not in gene_boundaries:
                raise ValueError(f"Gene {real_gene} not found in reference file.")

            b = gene_boundaries[real_gene.upper()]

        alignment = SeqIO.parse(input.alignment, "fasta")
        with open(output.alignment, "w") as oh:
            for record in alignment:
                sequence = Seq(record.seq)
                gene_keep = sequence[b[0]:b[1]]
                if set(gene_keep) in [{"N"}, {"-"}, set()]:
                    continue  # Skip sequences that are entirely masked
                sequence = len(sequence) * "-"
                sequence = sequence[:b[0]] + gene_keep + sequence[b[1]:]
                record.seq = Seq(sequence)
                SeqIO.write(record, oh, "fasta")

##############################
# Tree building
###############################
rule tree:
    message:
        """
        Segment: {wildcards.seg}
        Creating a maximum likelihood tree
        """
    input:
        # alignment = rules.align.output.alignment
        alignment = rules.sub_alignments.output.alignment
    benchmark:
        "benchmark/tree.{seg}{gene}{protein}.log"
    output:
        # tree = "{seg}/results/tree_raw.nwk"
        tree = "{seg}/results/tree_raw{gene}{protein}.nwk"
    threads: 9
    benchmark:
        "logs/tree.{seg}{gene}{protein}.log"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

##############################
# Refine tree &
# Ancestral sequence reconstruction
# & Translation
###############################
rule refine:
    message:
        """
        Segment: {wildcards.seg}
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
          - rooting with {params.rooting}
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.sub_alignments.output.alignment,
        # alignment = rules.align.output.alignment,
        metadata =  rules.add_metadata.output.metadata,
    benchmark:
        "benchmark/refine.{seg}{gene}{protein}.log"
    output:
        # tree = "{seg}/results/tree.nwk",
        # node_data = "{seg}/results/branch_lengths.json"
        tree = "{seg}/results/tree{gene}{protein}.nwk",
        node_data = "{seg}/results/branch_lengths{gene}{protein}.json",
    benchmark:
        "logs/refine.{seg}{gene}{protein}.log"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 10, # originally 3; set to 6 if you want more control over outliers
        strain_id_field = config["id_field"],
        clock_rate = 0.004, # remove for estimation
        clock_std_dev = 0.0015,
        rooting = lambda wildcards: "--root mid_point" if wildcards.seg == "P1" else "--root mid_point"
        # rooting = F"--root {REFERENCE}",
        # rooting = "--root AB575916 AB575939", ## B and C sequence
        # rooting = lambda wildcards: (
        #     "--root DQ341364 KF501389" if (wildcards.seg == "whole_genome" and not wildcards.gene)
        #     else "--root JN204010 DQ341364" if wildcards.seg == "vp1"
        #     else ""
        # )

    log:
        "logs/refine.{seg}{gene}{protein}.log" # number of dropped sequences
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
            --stochastic-resolve\
            --coalescent {params.coalescent} \
            --date-confidence \
            --clock-rate {params.clock_rate}\
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            {params.rooting} \
            2>&1 | tee {log}

        echo "Refined tree saved to {output.tree}"
        """
       

rule ancestral:
    message: 
        """
        Reconstructing ancestral sequences and mutations

        Segment: {wildcards.seg}
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.sub_alignments.output.alignment,
        # alignment = rules.align.output.alignment,
        annotation = rules.extract.output.extracted_genbank,
    output:
        node_data = "{seg}/results/muts{gene}{protein}.json",
    params:
        inference = "joint",
        genes = lambda wildcards: (
            CODING_GENES if wildcards.seg == "whole_genome"
            else PROTEIN1 if wildcards.seg == "P1"
            else "VP1" if wildcards.seg == "vp1"
            else []
        ),        
        translation_template= r"{seg}/results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"{seg}/results/translations/cds_%GENE.ancestral.fasta",
        root = "{seg}/results/ancestral_sequences.fasta",
    log:
        "logs/ancestral.{seg}{gene}{protein}.log"
    shell:
        """
            augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template} \
            --output-sequences {params.root} \
            > {log} 2>&1
        """
        # --keep-ambiguous\ #do not infer nucleotides at ambiguous (N) sites on tip sequences (leave as N).
        # --root-sequence {input.annotation} \  -> assigns mutations to the root relative to the reference, not wanted here
 
rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_metadata.output.metadata
    output:
        # node_data = "{seg}/results/traits.json"
        node_data = "{seg}/results/traits{gene}{protein}.json",
    params:
        traits = "country",
        strain_id_field= config["id_field"]
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

##############################
# Clade assignment
###############################

rule clades: 
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        # clade_data = "{seg}/results/clades.json"
        clade_data = "{seg}/results/clades{gene}{protein}.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule clade_published:
    message: "Assigning clades from publications for {wildcards.seg}"
    input:
        metadata     = rules.add_metadata.output.metadata,
        subgenotypes = files.VP1_CLADES,
        rivm_data    = files.RIVM_CLADES,
        c1           = "data/C1_colors_P1.csv"
    params:
        strain_id_field = config["id_field"],
        blast           = "{seg}/results/blast_{seg}_length.csv",

    output:
        final_metadata = "{seg}/results/final_metadata.tsv"
    run:
        import pandas as pd

        sid = params.strain_id_field

        meta             = pd.read_csv(input.metadata,     sep='\t', index_col=False)
        rivm_subtypes    = pd.read_csv(input.rivm_data,    index_col=False)
        vp1_subgenotypes = pd.read_csv(input.subgenotypes, sep='\t')
        blast_length     = pd.read_csv(params.blast, names=["accession", "length"])
        C1_colors        = pd.read_csv(input.c1)

        meta          = meta.drop_duplicates(subset=sid)
        rivm_subtypes = rivm_subtypes.drop_duplicates(subset="name")

        meta = pd.merge(meta, vp1_subgenotypes, on=sid, how="left")

        mask = (
            (rivm_subtypes["VP1 subgenogroup"] != "Could not assign") &
            (rivm_subtypes["type"] == "EV-A71")
        )
        rivm_subtypes = rivm_subtypes.loc[mask, ["name", "VP1 subgenogroup"]].copy()
        rivm_subtypes.rename(
            columns={"name": sid, "VP1 subgenogroup": "RIVM_subgenogroup"},
            inplace=True
        )
        meta = pd.merge(meta, rivm_subtypes, on=sid, how="left")
        meta["RIVM_subgenogroup"] = meta["RIVM_subgenogroup"].str.strip()

        meta = pd.merge(meta, blast_length, left_on=sid, right_on="accession", how="left")
        if sid != "accession" and "accession" in meta.columns:
            meta.drop(columns=["accession"], inplace=True)

        meta["length_blast"] = meta["length"].astype(str)
        meta.drop(columns=["length"], inplace=True)

        meta["url"] = "https://www.ncbi.nlm.nih.gov/nuccore/" + meta[sid]

        meta.rename(columns={'has_age':'Age available'}, inplace=True)
        meta.rename(columns={'has_diagnosis':'Diagnosis available'}, inplace=True)

        meta = pd.merge(meta, C1_colors, left_on=sid, right_on="accession", how="left")

        meta.to_csv(output.final_metadata, sep='\t', index=False)

##############################
# Assign epitopes to the tree colors
###############################
rule epitopes:
    input:
        anc_seqs = rules.ancestral.output.node_data, #"results/muts_vp1.json",
        tree = rules.refine.output.tree #"results/tree_vp1.nwk"
    output:
        node_data = "{seg}/results/epitopes{gene}{protein}.json"
    params:
        translation = "vp1/results/translations/cds_VP1.ancestral.fasta",
        epitopes = { #it is minus one!!!!
        'BC':     list(range(95, 107)),         # Huang et al., 2015; Foo et al., 2008; structural mapping of neutralizing antibodies
        'DE':     list(range(142, 152)),        # Liu et al., 2011; Zaini et al., 2012.
        'EF':     list(range(165, 173)),        # Lyu et al., 2014; Wang et al., 2010.
        'CTERM':  list(range(281, 291)),        # Chang et al., 2012 (monoclonal antibody studies), structural models.
        'GH':     list(range(209, 224)),        # mutations at S215, K218 have been noted to impact neutralization. Often used in vaccine design (e.g., in VLPs or epitope grafting studies).      
        'Esc_CHN':[283, 293]},                  # Escape Mutations in C-Terminal in Chinese Samples
        min_count = 6 # number of sequences?
    run:
        import json
        from collections import defaultdict
        from Bio import Phylo, SeqIO

        manyXList = ["XXXXXXXXXXXX", "KEXXXXXXXXXX", "KERANXXXXXXX", "KERXXXXXXXXX", "KERAXXXXXXXX"]
        valid_esc_chn = {"SA", "TA", "TX", "TS", "SS", "XX"}  # Set of valid values for Esc_CHN

        # Read translation files
        vp1_anc = SeqIO.to_dict(SeqIO.parse(params.translation, "fasta"))

        T = Phylo.read(input.tree, 'newick')
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        nodes = {}
        epitope_counts = {epi: defaultdict(int) for epi in params.epitopes}

        for node in T.find_clades(order='preorder'):
            n = node.name
            aa = vp1_anc[n].seq
            nodes[n] = {}
            for epi,pos in params.epitopes.items():
                pos = [p - 1 for p in pos]  # Convert to 0-based indexing
                nodes[n][epi] = "".join([aa[p] for p in pos])
                if epi == 'CTERM':
                    if nodes[n]['CTERM'] in manyXList:
                        nodes[n]['CTERM'] = "many x"
                    elif 'X' in nodes[n]['CTERM']:
                        nodes[n]['CTERM'] = nodes[node.parent.name]['CTERM']
                if epi == 'Esc_CHN':
                    if nodes[n]['Esc_CHN'] not in valid_esc_chn:
                        nodes[n]['Esc_CHN'] = "other"
                if not n.startswith('NODE_'):
                    epitope_counts[epi][nodes[n][epi]] += 1

        for node in nodes:
            for epi,seq in nodes[node].items():
                min_count2 = params.min_count if epi != "CTERM" else 6
                if epi == "CTERM" and seq in manyXList:
                    nodes[node][epi]='many X'
                elif epitope_counts[epi][seq]<min_count2:#params.min_count:
                    nodes[node][epi]='other'

        with open(output.node_data, 'w') as fh:
            json.dump({"epitopes": params.epitopes, "nodes":nodes}, fh)

##############################
# Assign key mutations to tree nodes for coloring in Nextstrain
# Output: node_data JSON for use with augur export --node-data
##############################

rule key_muts:
    input:
        anc_seqs = rules.ancestral.output.node_data,
        tree = rules.refine.output.tree,
        gff = files.gff_reference
    output:
        node_data = "{seg}/results/key_muts{gene}{protein}.json"
    params:
        muts = {
            "C1": {
                "VP1": ["16M","19V", "119L","262V", "262I", "289T"],
                "VP2": ["170V", "170I", "195M", "219V"],
                "VP3": ["93D", "232S"],
                "VP4": ["54T"],
            },
            "C2": {
                "VP1": ["22R", "31D", "145Q"],
                "VP2": ["144T", "170V", "224Y"],
                "VP3": ["62D", "93S", "155I"],
                "VP4": ["7A"],
            },
        },
        min_count = 1
    run:
        import json
        from collections import defaultdict
        from Bio import Phylo
        from augur.translate import safe_translate

        def mut_to_pos(mut_str: str) -> int:
            digits = "".join(ch for ch in mut_str if ch.isdigit())
            if not digits:
                raise ValueError(f"No position found in mutation string: {mut_str!r}")
            return int(digits)

        def parse_gff_cds_coords(gff_path: str) -> dict[str, tuple[int, int]]:
            coords = {}
            with open(gff_path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    seqid, source, ftype, start, end, score, strand, phase, attrs = line.split("\t")
                    if ftype != "CDS":
                        continue

                    # attributes are semicolon-separated key=value
                    attr_map = {}
                    for item in attrs.split(";"):
                        if "=" in item:
                            k, v = item.split("=", 1)
                            attr_map[k] = v

                    name = attr_map.get("Name")
                    if not name:
                        continue
                    # We only want the protein sub-CDSs, not the full P1 CDS
                    if name not in {"VP1", "VP2", "VP3", "VP4"}:
                        continue

                    coords[name] = (int(start), int(end))

            missing = {"VP1", "VP2", "VP3", "VP4"} - set(coords)
            if missing:
                raise ValueError(f"GFF missing CDS coords for: {sorted(missing)}")
            return coords

        def protein_aa_to_p1_aa_pos0(protein: str, aa_pos_in_protein_1based: int, cds_coords: dict[str, tuple[int, int]]) -> int:
            start_nt_1based, end_nt_1based = cds_coords[protein]
            protein_len_aa = (end_nt_1based - start_nt_1based + 1) // 3
            if not (1 <= aa_pos_in_protein_1based <= protein_len_aa):
                raise ValueError(
                    f"{protein} AA position {aa_pos_in_protein_1based} out of range (1..{protein_len_aa}) "
                    f"based on GFF coords {start_nt_1based}..{end_nt_1based}"
                )

            p1_start_aa0 = (start_nt_1based - 1) // 3
            return p1_start_aa0 + (aa_pos_in_protein_1based)

        with open(input.anc_seqs) as fh:
            anc = json.load(fh)["nodes"]

        cds_coords = parse_gff_cds_coords(input.gff)

        T = Phylo.read(input.tree, 'newick')
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        pos0_by_feature = {}
        feature_order = []
        for mapping_name, proteins in params.muts.items():
            feature_order.append(mapping_name)

            seen = set()
            pos0_list = []

            # iterate proteins in a stable order (optional, but keeps output reproducible)
            for protein in sorted(proteins.keys()):
                for mut in proteins[protein]:
                    aa_pos_in_protein = mut_to_pos(mut)
                    pos0 = protein_aa_to_p1_aa_pos0(protein, aa_pos_in_protein, cds_coords)
                    if pos0 not in seen:
                        pos0_list.append(pos0)
                        seen.add(pos0)

            pos0_by_feature[mapping_name] = pos0_list
        ## flatten pos0_by_feature to two list of positions (C1 and C2) and check for duplicates

        nodes = {}
        epitope_counts = {feature: defaultdict(int) for feature in feature_order}

        for node in T.find_clades(order="preorder"):
            n = node.name
            aa = safe_translate(anc[n]["sequence"])
            nodes[n] = {}

            for feature in feature_order:  # feature is "C1" or "C2"
                pos0_list = pos0_by_feature[feature]
                pos0_list.sort()
                nodes[n][feature] = "".join(aa[p-1] for p in pos0_list)

                if not n.startswith("NODE_"):
                    epitope_counts[feature][nodes[n][feature]] += 1

        for n in nodes:
            for feature, seq in nodes[n].items():
                if epitope_counts[feature][seq] < params.min_count:
                    nodes[n][feature] = "other"

        with open(output.node_data, "w") as fh:
            json.dump({"epitopes": params.muts, "nodes": nodes}, fh)


#########################
#  EXPORT
#########################
rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.clade_published.output.final_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        # nt_muts = "{seg}/results/nt_muts.json",
        # aa_muts = "{seg}/results/aa_muts.json",
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        vaccine = "config/vaccine.json",
        auspice_config = files.auspice_config,
        config_dates = "config/date_bounds.json",
        muts = rules.ancestral.output.node_data,
        epis = lambda wildcards: "vp1/results/epitopes.json" if wildcards.seg == "vp1" else [], ## please run the epitopes function
        key_muts =  lambda wildcards: "P1/results/key_muts.json" if wildcards.seg == "P1" else [],
    params:
        strain_id_field= config["id_field"],
        muts_flag = lambda wildcards: f"{wildcards.seg}/results/muts.json" if wildcards.seg else "",
        europe_col = lambda wildcards: "config/EUROPE_colors.tsv" if wildcards.seg == "P1" else files.colors,
    benchmark:
        "benchmark/export.{seg}{gene}{protein}.log"
    output:
        auspice_json="auspice/enterovirus_A71_{seg}{gene}{protein}.json"
        # auspice_json = "auspice/enterovirus_A71_{seg}-accession.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {params.muts_flag} {input.clades} \
                {input.vaccine} {input.epis} {input.key_muts} \
            --colors {params.europe_col}\
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
        # {input.epis} 
        
rule rename_whole_genome:
    message: "Rename whole-genome built"
    input: 
        json="auspice/enterovirus_A71_whole_genome.json"
    output:
        json="auspice/enterovirus_A71_whole-genome.json" # easier view in auspice
    shell:
        """
        mv {input.json} {output.json}
        """


rule rename_genes:
    message: 
        "Rename and compress the single genome builts"
    input: 
        json="auspice/enterovirus_A71_whole_genome{gene}.json"

    output:
        json="auspice/enterovirus_A71_gene_{gene}.json" # easier view in auspice
    shell:
        """
        jq -c . {input.json} > {output.json}
        rm {input.json}
        """


rule rename_proteins:
    message: 
        "Rename the single genome builts"
    input: 
        json="auspice/enterovirus_A71_whole_genome{protein}.json"

    output:
        json="auspice/enterovirus_A71_protein_{protein}.json"
    shell:
        """
        mv {input.json} {output.json}
        """


rule clean:
    message: "Removing directories: {params}"
    params:
        "ingest/data/*.*",
        "*/results/*",
        "auspice/*.json",
        "temp/*",
        "logs/*",
        "benchmark/*",
        files.METADATA,
        files.SEQUENCES,
        "data/curated/*",
        "data/all_sequences.fasta",
        "data/all_metadata.tsv",
        "data/final_metadata.tsv",
        "logs/*"
    shell:
        "rm -rfv {params}"


rule upload: ## make sure you're logged in to Nextstrain
    message: "Uploading auspice JSONs to Nextstrain"
    input:
        jsons = ["auspice/enterovirus_A71_vp1.json", "auspice/enterovirus_A71_whole-genome.json"]
        # "auspice/enterovirus_A71_gene_-vp1.json", "auspice/enterovirus_A71_gene_-3D.json"]
    params:
        remote_group=REMOTE_GROUP,
        date=UPLOAD_DATE,
    shell:
        """
        nextstrain login
        nextstrain remote upload \
            nextstrain.org/groups/{params.remote_group}/ \
            {input.jsons}
        # nextstrain logout
        mkdir -p auspice/{params.date}
        cp {input.jsons} auspice/{params.date}/
        """
