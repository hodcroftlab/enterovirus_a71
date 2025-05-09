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

# Load config file
if not config:
    configfile: "config/config.yaml"

###############
wildcard_constraints:
    seg="vp1|whole_genome",
    gene="|-5utr|-vp4|-vp2|-vp3|-vp1|-2A|-2B|-2C|-3A|-3B|-3C|-3D|-3utr",
    protein="|-P1|-P2|-P3"
   
# Define segments to analyze
segments = ['vp1', 'whole-genome']
GENES=["-5utr","-vp4", "-vp2", "-vp3", "-vp1", "-2A", "-2B", "-2C", "-3A", "-3B", "-3C", "-3D","-3utr"]
PROT = ["-P1", "-P2", "-P3"]

# include or exclude ENPEN data
INCL_ENPEN = False

# Expand augur JSON paths
rule all:
    input:
        augur_jsons = expand("auspice/enterovirus_A71_{segs}.json", segs=segments),
        epitopes = "vp1/results/epitopes.json"

rule all_genes:
    input:
        seg_jsn = expand("auspice/enterovirus_A71_{segs}.json", segs=segments),
        augur_jsons = expand("auspice/enterovirus_A71_gene_{genes}.json", genes=GENES),
        epitopes = "vp1/results/epitopes.json"
        

rule all_proteins:
    input:
        augur_jsons = expand("auspice/enterovirus_A71_protein_{proteins}.json", proteins=PROT),
        seg_jsn = expand("auspice/enterovirus_A71_{segs}.json", segs=segments),
        epitopes = "vp1/results/epitopes.json"


# Rule to handle configuration files
rule files:
    input:
        sequence_length =   "{seg}",
        dropped_strains =   "config/dropped_strains.txt",
        incl_strains =      "config/kept_strains.txt",
        reference =         "{seg}/config/reference_sequence.gb",
        gff_reference =     "{seg}/config/annotation.gff3",
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv",
        meta=               "data/metadata.tsv",
        extended_metafile=  "data/meta_public.tsv",
        meta_collab =       "data/meta_ENPEN.tsv",
        last_updated_file = "data/date_last_updated.txt",
        local_accn_file =   "data/local_accn.txt"


files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

# rule fetch:
#     input:
#         dir = "ingest"
#     output:
#         sequences="data/sequences.fasta",
#         metadata=files.meta
#     params:
#         seq="ingest/data/sequences.fasta",
#         meta="ingest/data/metadata.tsv"
#     shell:
#         """
#         cd {input.dir} 
#         snakemake --cores 9 all
#         cd ../
#         cp -u {params.seq} {output.sequences}
#         cp -u {params.meta} {output.metadata}
#         """

##############################
# Update strain names
###############################

# rule update_strain_names:
#     message:
#         """
#         Updating strain name in metadata.
#         """
#     input:
#         file_in =  files.meta
#     params:
#         backup = "data/strain_names_previous_run.tsv"
#     output:
#         file_out = "data/updated_strain_names.tsv"
#     shell:
#         """
#         time bash scripts/update_strain.sh {input.file_in} {params.backup} {output.file_out}
#         cp {output.file_out} {params.backup}
#         """


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
        metadata=files.extended_metafile,  # Path to input metadata file
        meta_collab = files.meta_collab  # Data shared with us by collaborators
    params:
        strain_id_field=config["id_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
    output:
        metadata = "data/curated/meta_public.tsv",  # Final output file for publications metadata
        meta_collab="data/curated/meta_collab.tsv",  # Curated collaborator metadata
        meta="data/curated/all_meta.tsv"  # Final merged output file
    shell:
        """
        # Normalize strings for publication metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {input.metadata} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        
        # Normalize strings and format dates for collab metadata
        augur curate normalize-strings \
            --id-column {params.strain_id_field} \
            --metadata {input.meta_collab} \
        | augur curate format-dates \
            --date-fields {params.date_fields} \
            --no-mask-failure \
            --expected-date-formats {params.expected_date_formats} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.meta_collab}
        
        # Merge curated metadata
        augur merge --metadata metadata={output.metadata} meta_collab={output.meta_collab}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.meta}
        """


##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################

rule update_sequences:
    input:
        sequences = "data/sequences.fasta",
        metadata=files.meta,
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
# blast fasta files for vp1 
###############################

rule blast:
    input: 
        blast_db_file = "data/references/reference_vp1_blast.fasta",
        seqs_to_blast = rules.update_sequences.output.sequences
    output:
        blast_out = "vp1/temp/blast_out.csv"
    params:
        blast_db = "vp1/temp/blast_database"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out {output.blast_out} -evalue 0.0005
        """

rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out, # output blast (for your protein)
        input_seqs = rules.update_sequences.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta"
        
    params:
        protein = [600,900], #TODO: min & max length for protein
        whole_genome = [6400,8000], #TODO: min & max length for whole genome
        range = "{seg}" # this is determining the path it takes in blast_sort (protein-specific or whole genome)
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --protein_length {params.protein}  --whole_genome_length {params.whole_genome} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --range {params.range}
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
        metadata=files.meta,
        new_data=rules.curate.output.meta,
        regions=ancient(files.regions),
        renamed_strains="data/updated_strain_names.tsv"
    params:
        strain_id_field=config["id_field"],
        last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
        C1like_accn = "data/list_C1like_accn.txt",
    output:
        metadata="data/all_metadata.tsv"
    shell:
        """
        python scripts/add_metadata.py \
            --input {input.metadata} \
            --add {input.new_data} \
            --rename {input.renamed_strains} \
            --local {params.local_accn} \
            --C1like {params.C1like_accn} \
            --update {params.last_updated}\
            --regions {input.regions} \
            --id {params.strain_id_field} \
            --output {output.metadata}
        
        if [ -d "./temp/" ]; then
        rm -r ./temp/
        fi
        """

##############################
# Create an index of sequence composition for filtering & filter
###############################
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
        Log file in {log.log}
        """
    input:
        sequences = rules.blast_sort.output.sequences, ## x had no sequence data -> dropped because they don't meet the min sequence length in blast_sort
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.add_metadata.output.metadata,
        exclude = files.dropped_strains,
        include = files.incl_strains,
    output:
        sequences = "{seg}/results/filtered.fasta",
        reason ="{seg}/results/reasons.tsv",
    log:
        log = "logs/filter.{seg}.log"
    benchmark:
        "benchmark/filter.{seg}.log"

    params:
        group_by = "country year",
        sequences_per_group = 5000, # 2000 originally
        strain_id_field=config["id_field"],
        min_date = 1960,  # BrCr was collected in 1970
        exclude_enpen = "--exclude-where ENPEN=True" if not INCL_ENPEN else "" # INCL_ENPEN was defined on line 32
    shell:
        """
        (augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            {params.exclude_enpen} \
            --include {input.include} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --output {output.sequences}\
            --output-log {output.reason}) > {log.log} 2>&1
        """

##############################
# Reference sequence &
# Alignment
###############################
rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format
        """
    input:
        reference = files.reference

    output:
        reference = "{seg}/results/reference_sequence.fasta"
    run:
        from Bio import SeqIO 
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextclade.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "{seg}/results/aligned.fasta"
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
        --output-fasta {output.alignment} 
        """

# potentially add one-by-one genes
# use wildcards
rule sub_alignments:
    input:
        alignment=rules.align.output.alignment,
        reference=files.reference
    output:
        # alignment = "{seg}/results/aligned.fasta"
        alignment = "{seg}/results/aligned{gene}{protein}.fasta"
    benchmark:
        "benchmark/sub_alignments.{seg}{gene}{protein}.log"
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq

        if wildcards.protein:
            real_gene = wildcards.protein.replace("-", "", 1)
            boundaries = {
                'P1':(744,3329),  'P2':(3330,5063),
                'P3':(5064,7322)}
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
                if set(gene_keep) == {"N"} or len(gene_keep) == 0 or set(gene_keep) == {"-"}:
                    continue  # Skip sequences that are entirely masked
                sequence = len(sequence) * "N"
                sequence = sequence[:b[0]] + gene_keep + sequence[b[1]:]
                record.seq = Seq(sequence)
                SeqIO.write(record, oh, "fasta")

##############################
# Tree building
###############################
rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        # alignment = rules.fix_align_codon.output.alignment
        alignment = rules.sub_alignments.output.alignment
    benchmark:
        "benchmark/tree.{seg}{gene}{protein}.log"
    output:
        # tree = "{seg}/results/tree_raw.nwk"
        tree = "{seg}/results/tree_raw{gene}{protein}.nwk"
    threads: workflow.cores
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
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
        metadata =  rules.add_metadata.output.metadata,
        reference = rules.reference_gb_to_fasta.output.reference
    benchmark:
        "benchmark/refine.{seg}{gene}{protein}.log"
    output:
        # tree = "{seg}/results/tree.nwk",
        # node_data = "{seg}/results/branch_lengths.json"
        tree = "{seg}/results/tree{gene}{protein}.nwk",
        node_data = "{seg}/results/branch_lengths{gene}{protein}.json",
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 3, # originally 3; set to 6 if you want more control over outliers
        strain_id_field = config["id_field"],
        clock_rate = 0.004, # remove for estimation
        clock_std_dev = 0.0015,
        rooting = lambda wildcards: (
            "--root DQ341364 KF501389" if (wildcards.seg == "whole_genome" and not wildcards.gene)
            else "--root JN204010 DQ341364" if wildcards.seg == "vp1"
            else ""
        )

    log:
        reasons_refine = "logs/refine.{seg}{gene}{protein}.log" # number of dropped sequences
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
            --clock-rate {params.clock_rate}\
            --clock-std-dev {params.clock_std_dev} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            {params.rooting} >> {log.reasons_refine} 2>&1
        """
        #         dropped=$(comm -23 \
        #     <(grep -oE '[^\(\),:]+' {input.tree} | sort -u) \
        #     <(grep -oE '[^\(\),:]+' {output.tree} | sort -u)
        # )
        # dropped_count=$(echo "$dropped" | wc -l)

        # echo "Dropped sequences due to clock filter: $dropped_count" >> {log.reasons_refine}
        # echo "Dropped tip labels:" >> {log.reasons_refine}
        # echo "$dropped" >> {log.reasons_refine}
        # echo -e "\\n Excluded tips count:" >> {log.reasons_refine}
        # grep "pruning leaf" {log.reasons_refine} | wc -l >> {log.reasons_refine}
        

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.sub_alignments.output.alignment
    output:
        # node_data = "{seg}/results/nt_muts.json"
        node_data = "{seg}/results/nt_muts{gene}{protein}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous\
            --inference {params.inference}
        """
 
rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "{seg}/results/aa_muts{gene}{protein}.json"
        # node_data = "{seg}/results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

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
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        # clade_data = "{seg}/results/clades.json"
        clade_data = "{seg}/results/clades{gene}{protein}.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule clade_published:
    message: "Assigning clades from publications"
    input:
        metadata = rules.add_metadata.output.metadata,
        subgenotypes = "data/clades_vp1.tsv",
        rivm_data = "data/rivm/subgenotypes_rivm.csv",
        alignment= "vp1/results/aligned.fasta"
    params:
        strain_id_field= config["id_field"],
        rerun=True
    output:
        final_metadata = "data/final_metadata.tsv"
    shell:
        """
        python scripts/published_clades.py --input {input.metadata} --rivm {input.rivm_data}\
        --sgt {input.subgenotypes} --alignment {input.alignment} --id {params.strain_id_field} --output {output.final_metadata}
        """

##############################
# Assign epitopes to the tree colors
###############################
rule epitopes:
    input:
        anc_seqs = rules.ancestral.output.node_data, #"results/nt_muts_vp1.json",
        tree = rules.refine.output.tree #"results/tree_vp1.nwk"
    output:
        node_data = "{seg}/results/epitopes{gene}{protein}.json"
    params:
        epitopes = {
            # 'BC':[94, 96, 97, 98, 99, 100, 101, 102, 103, 104], 
            # 'EF':{165, 167, 168},
            # 'DE':[141,142,143,144,145,146,147,148,149],
            # 'CTERM':[272, 273, 280, 281, 282, 283, 292,293], 
            'Esc_CHN':[282, 292]}, # it is minus one!!!!
        min_count = 6 # number of sequences?
    run:
        import json
        from collections import defaultdict
        from augur.translate import safe_translate
        from Bio import Phylo

        manyXList = ["XXXXXXXXXXXX", "KEXXXXXXXXXX", "KERANXXXXXXX", "KERXXXXXXXXX", "KERAXXXXXXXX"]
        valid_esc_chn = {"SA", "TA", "TX", "TS", "SS", "XX"}  # Set of valid values for Esc_CHN
        with open(input.anc_seqs) as fh:
            anc = json.load(fh)["nodes"]

        T = Phylo.read(input.tree, 'newick')
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        nodes = {}
        epitope_counts = {epi: defaultdict(int) for epi in params.epitopes}

        for node in T.find_clades(order='preorder'):
            n = node.name
            aa = safe_translate(anc[n]["sequence"])

            nodes[n] = {}
            for epi,pos in params.epitopes.items():
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
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        vaccine = "config/vaccine.json",
        auspice_config = files.auspice_config
    params:
        strain_id_field= config["id_field"],
        epis = lambda wildcards: "vp1/results/epitopes.json" if wildcards.seg == "vp1" else "", ## please run the epitopes function
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
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
                {input.aa_muts} {input.clades} {input.vaccine} {params.epis} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
        # {input.epis} 
        
rule clean:
    message: "Removing directories: {params}"
    params:
        "results/*",
        # "auspice/*",
        "temp/*"
    shell:
        "rm -rfv {params}"


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
        "Rename the single genome builts"
    input: 
        json="auspice/enterovirus_A71_whole_genome{gene}.json"

    output:
        json="auspice/enterovirus_A71_gene_{gene}.json" # easier view in auspice
    shell:
        """
        mv {input.json} {output.json}
        """


rule rename_proteins:
    message: 
        "Rename the single genome builts"
    input: 
        json="auspice/enterovirus_A71_whole-genome{protein}.json"

    output:
        json="auspice/enterovirus_A71_protein_{protein}.json"
    shell:
        """
        mv {input.json} {output.json}
        """