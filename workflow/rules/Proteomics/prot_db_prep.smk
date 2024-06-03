### Proteomics search database preparation ###
import glob

rule prep_genome_updater:
    output:
        os.path.join(ENVDIR, 'genome_updater.ready')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(4)
    conda:
        os.path.join(ENVDIR, "IMP_protein_db.yaml")
    log:
         os.path.join(ENVDIR, "logs/prep_genome_updater.log")
    script:
        os.path.join(SRCDIR, "prep_genome_updater.py")

if config['proteomics']['add_kraken_taxa_proteomes']:
    # Get genome and proteome of top representatives of all species and strains (tax level S and S1) identified confidently by kraken
    checkpoint download_kraken_tax_proteomes:
        input:
            genome_updater_check = os.path.join(ENVDIR, 'genome_updater.ready'),
            kraken_reports = expand("Analysis/taxonomy/kraken/{type}.kraken.report", type=TYPES)
        output:
            proteomics_dir = directory("Analysis/annotation/proteomics/taxa_proteins"),
            download_dir = directory("Analysis/annotation/proteomics/taxa_proteins/download"),
            download_summary = "Analysis/annotation/proteomics/taxa_proteins/assembly_summary.txt"
        resources:
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            getThreads(128)
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/download_kraken_tax_proteomes.log"
        shell:
            """
            echo "Gathering all kraken taxa identified ..."
            all_taxa=$(grep -hP "\tS\t|\tS1\t" {input.kraken_reports} | cut -f 7 | sort | uniq | tr '\n' ',' | sed 's/,$//')
            echo "Found $(echo $all_taxa | tr ',' '\n' | wc -l) taxa." >> {log}

            echo "Downloading genomes and proteomes for all taxa ..." >> {log}
            genome_updater.sh -M "ncbi" \
                              -A 1 \
                              -t {threads} \
                              -L curl \
                              -d "refseq,genbank" \
                              -f genomic.fna.gz,protein.faa.gz \
                              -b download \
                              -T $all_taxa \
                              -o {output.proteomics_dir} >> {log} 2>&1

            echo "Done." >> {log}
            """

    # Check if proteomes for taxa are available and valid fastas, if not predict proteins (unless eukaryote)
    rule process_kraken_taxa_proteomes:
        input:
            # `get_kraken_tax_genomes` returns {download_kraken_tax_proteomes.output.proteomics_dir}/download/files/{genome_id}_genomic.fna.gz
            # kraken_taxa_genome = get_kraken_tax_genomes,
            kraken_taxa_genome = "Analysis/annotation/proteomics/taxa_proteins/download/files/{genome_id}_genomic.fna.gz",
            download_dir = "Analysis/annotation/proteomics/taxa_proteins/download",
            assembly_summary = "Analysis/annotation/proteomics/taxa_proteins/assembly_summary.txt"
        output:
            process_check = "status/proteomics/{genome_id}.processed"
        resources:
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            genome_id = lambda wildcards: wildcards.genome_id,
            proteome_fasta = lambda wildcards: os.path.join("Analysis/annotation/proteomics/taxa_proteins/download/files", f"{wildcards.genome_id}_protein.faa.gz"),
            unsupported_taxa = 'plant|vertebrate_mammalian|fungi|invertebrate|protozoa|vertebrate_other'
        threads:
            getThreads(128)
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/prep_kraken_genomes/process_{genome_id}_proteome.log"
        shell:
            """
            echo "{params.genome_id}: Processing genome/proteome ..." >> {log}

            # Check if the proteome file exists and is valid. if `seqkit stats` doesnt fail we dont need to process further
            echo "{params.genome_id}: Checking proteome fasta ..." >> {log}
            seqkit_check=$(seqkit stats --seq-type protein "{params.proteome_fasta}" 2>&1)
            if echo "$seqkit_check" | grep -qv 'invalid FASTA/Q format'; then
                echo "{params.genome_id}: Valid proteome fasta. No gene prediction needed. Exiting." >> {log}
                touch {output.process_check}
                exit 0
            fi

            # Exit if the proteome file is not a valid fasta file and the genome is an euakryote
            echo "{params.genome_id}: Checking if genome is eukaryote ..." >> {log}
            if [[ $(grep {params.genome_id} {input.assembly_summary} | grep -q -E "{params.unsupported_taxa}" ; echo $?) -eq 0 ]]; then
                echo "{params.genome_id}: Eukaryote without valid proteome fasta. Deleting proteome file (if it exists) and exiting." >> {log}
                rm -rf {params.proteome_fasta} >> {log} 2>&1
                touch {output.process_check}
                exit 0
            fi

            # If the proteome fasta is not valid or not present, run gene prediction
            echo "{params.genome_id}: Running Pyrodigal ..." >> {log}
            python3 {SRCDIR}/predict_proteins.py --input {input.kraken_taxa_genome} \
                                                 --output {params.proteome_fasta} \
                                                 --cores {threads} \
                                                 --write_gff >> {log} 2>&1
            touch {output.process_check}
            echo "{params.genome_id}:Done." >> {log}
            """

    # Add all kraken taxa proteomes to the bakta metaproteome
    rule add_kraken_taxa_proteomes:
        input: 
            bakta_fasta = "Analysis/annotation/bakta/annotation.faa",
            download_dir = "Analysis/annotation/proteomics/taxa_proteins/download",
            kraken_proteomes_processed = gather_kraken_tax_proteome_checks
        output:
            combined_fasta = temp("Analysis/annotation/proteins_kraken_taxa.faa")
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/add_kraken_taxa_proteomes.log"
        shell:
            """
            echo "Adding kraken taxa proteomes to the bakta protein fasta ..." >> {log}
            # Find of remaining proteome fastas
            kraken_proteomes=$(find {input.download_dir}/files -name "*_protein.faa.gz")
            cat {input.bakta_fasta} <(zcat $kraken_proteomes) > {output.combined_fasta}
            echo "Done." >> {log}
            """

    # Set kraken taxa proteome combined fasta as input for further processing
    PROTEIN_FASTA_RAW = "Analysis/annotation/proteins_kraken_taxa.faa"
else:
    # Set bakta metaproteome as input for further processing
    PROTEIN_FASTA_RAW = "Analysis/annotation/bakta/annotation.faa"


### Add extra proteomes, if specified ###
if config['proteomics']['extra_proteomes'].split():
    # Merge host proteome with metaproteome
    rule add_extra_proteomes:
        input:
            PROTEIN_FASTA_RAW
        output:
            temp("Analysis/annotation/proteins_host.faa")
        resources: 
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            extra_proteomes = config['proteomics']['extra_proteomes'].split()
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/add_extra_proteomes.log"
        message:
            "Adding extra proteome(s) to the metaproteomics proteome ..."
        shell:
            """
            echo "Merging extra proteome(s) with metaproteomics proteome ..." >> {log}
            cat {input} {params.extra_proteomes} > {output} 2>> {log}
            echo "Done." >> {log}
            """
    
    # Set extra proteome merged fasta as input for further processing
    PROTEIN_FASTA_PRE_DEDUP = "Analysis/annotation/proteins_host.faa"
else:
    # Set bakta metaproteome as input for further processing
    PROTEIN_FASTA_PRE_DEDUP = PROTEIN_FASTA_RAW

### Remove duplicate proteins ###
rule process_protein_fasta:
    input:
        protein_fasta = PROTEIN_FASTA_PRE_DEDUP
    output:
        dedup_fasta = temp("Analysis/annotation/proteins_dedup.faa.tmp"),
        duplicate_list = temp("Analysis/annotation/protein_duplicates.txt.tmp")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        max_len = config['proteomics']['max_len'],
        min_len = config['proteomics']['min_len'],
        sample = config['sample']
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_protein_db.yaml")
    log:
        "logs/analysis/remove_duplicate_proteins.log"
    shell:
        """
        echo "Filtering proteins by size, reformatting headers and removing duplicate proteins  ..." >> {log}
        seqkit seq --max-len {params.max_len} \
                   --min-len {params.min_len} \
                   --line-width 0 {input.protein_fasta} 2>> {log} | \
            sed -e "s/>/>{params.sample}___/g" \
                -e "s/ /__/g" 2>> {log} | \
            seqkit rmdup --by-seq \
                         --dup-num-file {output.duplicate_list} \
                         > {output.dedup_fasta} 2>> {log}
        echo "Done." >> {log}
        """

# Rename proteins to simple numerical IDs and create a mapping file to original IDs
rule rename_protein_fasta:
    input:
        dedup_fasta = "Analysis/annotation/proteins_dedup.faa.tmp",
        duplicate_list = "Analysis/annotation/protein_duplicates.txt.tmp"
    output:
        renamed_fasta = temp("Analysis/annotation/proteins_dedup.faa"),
        id2duplicates = "Analysis/annotation/protein_id2duplicates.tsv"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        sample = config['sample']
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_protein_db.yaml")
    log:
        "logs/analysis/rename_proteome_duplicates.log"
    shell:
        """
        echo "Renaming duplicate proteins ..." >> {log}
        python3 {SRCDIR}/rename_sample_rmdup_fasta.py {params.sample} \
                                                      {input.dedup_fasta} \
                                                      {input.duplicate_list} \
                                                      {output.renamed_fasta} \
                                                      {output.id2duplicates} >> {log} 2>&1
        echo "Done." >> {log}
        """

# Add CRAP proteins to the protein fasta
rule add_crap:
    input:
        renamed_fasta = "Analysis/annotation/proteins_dedup.faa"
    output:
        final_fasta = "Analysis/annotation/proteomics.final.faa.gz"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        crap_fasta = config['proteomics']['crap_db']
    threads:
        getThreads(12)
    log:
        "logs/analysis/add_crap_to_fasta.log"
    shell:
        """
        echo "Adding CRAP proteins to the protein fasta ..." >> {log}
        cat {input.renamed_fasta} \
            {params.crap_fasta} 2>> {log} \
          | sed -e "s/ /__/g" 2>> {log} \
          | pigz -p {threads} -c > {output.final_fasta} 2>> {log}
        echo "Done." >> {log}
        """
