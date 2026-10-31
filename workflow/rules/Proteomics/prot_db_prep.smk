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
            download_dir = directory("Analysis/annotation/proteomics/taxa_proteins/download"),
            download_files_dir = temp(directory("Analysis/annotation/proteomics/taxa_proteins/download/files")),
            download_summary = "Analysis/annotation/proteomics/taxa_proteins/assembly_summary.txt",
            proteomics_dir = directory("Analysis/annotation/proteomics/taxa_proteins"),
        params:
            min_reads = config['proteomics']['min_read_support']
        resources:
            runtime = SHORTISH_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            getThreads(64)
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/download_kraken_tax_proteomes.log"
        shadow: "copy-minimal"
        shell:
            """
            echo "Gathering all kraken taxa identified ..." >> {log}
            all_taxa=$(awk -F $'\t' -v min_reads={params.min_reads} '$6 ~ /^S[1]?$/ && $2 >= min_reads {{ print $7 }}' {input.kraken_reports} | sort -u | tr '\\n' ',' | sed 's/,$//')
            echo "Found $(echo $all_taxa | tr ',' '\n' | wc -l) taxa." >> {log}

            echo "Downloading genomes and proteomes for all taxa ..." >> {log}
            genome_updater.sh -M "ncbi" \
                            -A 1 \
                            -t {threads} \
                            -L curl \
                            -d "refseq,genbank" \
                            -f genomic.fna.gz,protein.faa.gz,genomic.gff.gz  \
                            -b download \
                            -T $all_taxa \
                            -o {output.proteomics_dir} >> {log} 2>&1

            echo "Done." >> {log}
            """

    # Check if proteomes for taxa are available and valid fastas, if not predict proteins (unless eukaryote)
    rule process_kraken_taxa_proteomes:
        input:
            kraken_taxa_genome = "Analysis/annotation/proteomics/taxa_proteins/download/files/{genome_id}_genomic.fna.gz",
            download_dir = "Analysis/annotation/proteomics/taxa_proteins/download",
            download_files_dir = "Analysis/annotation/proteomics/taxa_proteins/download/files",
            assembly_summary = "Analysis/annotation/proteomics/taxa_proteins/assembly_summary.txt"
        output:
            process_check = "status/proteomics/{genome_id}.processed"
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            genome_id = lambda wildcards: wildcards.genome_id,
            genome_accession = lambda wildcards: '_'.join(wildcards.genome_id.split('_')[:2]),
            taxon_base_path = "Analysis/annotation/proteomics/taxa_proteins/download/files",
            unsupported_taxa = 'plant|vertebrate_mammalian|fungi|invertebrate|protozoa|vertebrate_other'
        threads:
            getThreads(16)
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/prep_kraken_genomes/process_{genome_id}_proteome.log"
        script:
            os.path.join(SRCDIR, "process_kraken_detected_genome.py")

    # rule process_kraken_taxa_proteomes:
    #     input:
    #         # `get_kraken_tax_genomes` returns {download_kraken_tax_proteomes.output.proteomics_dir}/download/files/{genome_id}_genomic.fna.gz
    #         # kraken_taxa_genome = get_kraken_tax_genomes,
    #         kraken_taxa_genome = "Analysis/annotation/proteomics/taxa_proteins/download/files/{genome_id}_genomic.fna.gz",
    #         download_dir = "Analysis/annotation/proteomics/taxa_proteins/download",
    #         assembly_summary = "Analysis/annotation/proteomics/taxa_proteins/assembly_summary.txt"
    #     output:
    #         process_check = "status/proteomics/{genome_id}.processed"
    #     resources:
    #         runtime = SHORTISH_RUNTIME,
    #         mem_mb = MAIN_MEMCORE_MB
    #     params:
    #         genome_id = lambda wildcards: wildcards.genome_id,
    #         genome_accession = lambda wildcards: '_'.join(wildcards.genome_id.split('_')[:2]),
    #         taxon_base_path = "Analysis/annotation/proteomics/taxa_proteins/download/files",
    #         proteome_fasta = lambda wildcards: os.path.join("Analysis/annotation/proteomics/taxa_proteins/download/files", f"{wildcards.genome_id}_protein.faa.gz"),
    #         unsupported_taxa = 'plant|vertebrate_mammalian|fungi|invertebrate|protozoa|vertebrate_other'
    #     threads:
    #         getThreads(16)
    #     conda:
    #         os.path.join(ENVDIR, "IMP_protein_db.yaml")
    #     log:
    #         "logs/analysis/prep_kraken_genomes/process_{genome_id}_proteome.log"
    #     shadow: "copy-minimal"
    #     shell:
    #         """
    #         echo "{params.genome_id}: Processing genome/proteome ..." >> {log}

    #         # Extract the URL from the assembly summary file
    #         url=$(awk -F'\t' -v acc="{params.genome_accession}" '($1==acc || $18==acc) {{print $20}}' {input.assembly_summary})
            
    #         if [ -z "$url" ]; then
    #             echo "{params.genome_id}: Unable to find URL in assembly summary. Exiting." >> {log}
    #             exit 1
    #         fi

    #         # Check if the downloaded archives are valid
    #         for archive in genomic.fna.gz protein.faa.gz genomic.gff.gz; do
    #             # Extract the filename from the URL
    #             filename=$(basename "$url")_$archive

    #             echo "{params.genome_id}: Checking $archive at {params.taxon_base_path}/$filename ..." >> {log}

    #             if ! gzip -t {params.taxon_base_path}/$filename 2>/dev/null; then
    #                 echo "{params.genome_id}: Invalid archive. Attempting to re-download." >> {log}
    #                 rm -rf {params.taxon_base_path}/$filename >> {log} 2>&1

    #                 # Try to download the file
    #                 set +e  # Temporarily disable exit on error
    #                 wget --tries=5 \
    #                     --wait=10 \
    #                     --random-wait \
    #                     --retry-connrefused \
    #                     --no-verbose \
    #                     --retry-on-http-error=404,500,503 \
    #                     -O {params.taxon_base_path}/$filename \
    #                     "$url/$filename" >> {log} 2>&1
    #                 set -e  # Re-enable exit on error

    #                 # Check if the download was successful
    #                 if ! gzip -t {params.taxon_base_path}/$filename 2>/dev/null; then
    #                     rm -rf {params.taxon_base_path}/$filename >> {log} 2>&1

    #                     echo "{params.genome_id}: Download failed. Attempting to get alternate accession from assembly summary." >> {log}

    #                     # Get alternate accession
    #                     alt_accession=$(awk -F'\t' -v acc="{params.genome_accession}" '$1==acc {{print $18}}' {input.assembly_summary})

    #                     # If alt_accession is empty, we cannot proceed
    #                     if [ -z "$alt_accession" ]; then
    #                         echo "{params.genome_id}: Alternate accession not found in assembly summary. Exiting." >> {log}
    #                         exit 1
    #                     fi

    #                     # Make alternate URL by swapping out the accession with the alternate accession
    #                     if [[ $alt_accession == GCF* ]]; then
    #                         alt_url=$(echo "$url" | sed -e "s/{params.genome_accession}/$alt_accession/" -e 's/GCA/GCF/')
    #                     elif [[ $alt_accession == GCA* ]]; then
    #                         alt_url=$(echo "$url" | sed -e "s/{params.genome_accession}/$alt_accession/" -e 's/GCF/GCA/')
    #                     else
    #                         echo "{params.genome_id}: Invalid alternate accession $alt_accession. Skipping to next archive." >> {log}
    #                         continue  # Skip to the next archive
    #                     fi

    #                     # Make alternate filename
    #                     alt_filename=$(basename "$alt_url")_$archive

    #                     echo "{params.genome_id}: Re-downloading with alternate URL $alt_url and accession $alt_accession ..." >> {log}
    #                     set +e  # Temporarily disable exit on error
    #                     wget --tries=5 \
    #                          --wait=10 \
    #                          --random-wait \
    #                          --retry-connrefused \
    #                          --no-verbose \
    #                          --retry-on-http-error=404,500,503 \
    #                          -O {params.taxon_base_path}/$filename \
    #                          "$alt_url/$alt_filename" >> {log} 2>&1
    #                     set -e  # Re-enable exit on error

    #                     # Check if the alternate download was successful
    #                     if gzip -t {params.taxon_base_path}/$filename 2>/dev/null; then
    #                         echo "{params.genome_id}: Successfully re-downloaded and validated archive {params.taxon_base_path}/$filename." >> {log}
    #                     else
    #                         echo "{params.genome_id}: Re-download with alternate accession failed or file {params.taxon_base_path}/$filename still invalid. Giving Up." >> {log}
    #                         rm -rf {params.taxon_base_path}/$filename >> {log} 2>&1
    #                     fi
    #                 else
    #                     echo "{params.genome_id}: Successfully re-downloaded and validated archive {params.taxon_base_path}/$filename." >> {log}
    #                 fi
    #             else
    #                 echo "{params.genome_id}: Archive {params.taxon_base_path}/$filename is valid." >> {log}
    #             fi
    #         done
            
    #         # Check if the proteome file exists and is valid. If `seqkit stats` doesn't fail, we don't need to process further
    #         echo "{params.genome_id}: Checking proteome fasta ..." >> {log}

    #         set +e  # Temporarily disable exit on error
    #         seqkit_check=$(set +o pipefail; seqkit stats --seq-type protein {params.proteome_fasta} 2>&1)
    #         seqkit_exit_code=$?
    #         set -e  # Re-enable exit on error

    #         if [ $seqkit_exit_code -eq 0 ] && ! echo "$seqkit_check" | grep -q 'invalid FASTA/Q format'; then
    #             echo "{params.genome_id}: Valid proteome fasta. No gene prediction needed. Exiting." >> {log}
    #             touch {output.process_check}
    #             exit 0
    #         else
    #             # Remove invalid proteome fasta
    #             echo "{params.genome_id}: Invalid proteome fasta. Deleting proteome file." >> {log}
    #             rm -rf {params.proteome_fasta} >> {log} 2>&1
    #         fi

    #         # Exit if the proteome file is not a valid fasta file and the genome is an eukaryote
    #         echo "{params.genome_id}: Checking if genome is eukaryote ..." >> {log}
    #         if grep -q {params.genome_id} {input.assembly_summary} | grep -q -E "{params.unsupported_taxa}"; then
    #             echo "{params.genome_id}: Eukaryote without valid proteome fasta. Exiting." >> {log}
    #             touch {output.process_check}
    #             exit 0
    #         fi

    #         # If the proteome fasta is not valid or not present, run gene prediction
    #         echo "{params.genome_id}: Running Pyrodigal ..." >> {log}
    #         python3 {SRCDIR}/predict_proteins.py --input {input.kraken_taxa_genome} \
    #                                             --output {params.proteome_fasta} \
    #                                             --cores {threads} \
    #                                             --write_gff >> {log} 2>&1
    #         touch {output.process_check}
    #         echo "{params.genome_id}: Done." >> {log}
    #         """

    # Add all kraken taxa proteomes to the bakta metaproteome
    rule add_kraken_taxa_proteomes:
        input: 
            bakta_fasta = "Analysis/annotation/bakta/annotation.faa",
            download_dir = "Analysis/annotation/proteomics/taxa_proteins/download",
            download_files_dir = "Analysis/annotation/proteomics/taxa_proteins/download/files",
            kraken_proteomes_processed = gather_kraken_tax_proteome_checks
        output:
            combined_fasta = temp("Analysis/annotation/proteins_kraken_taxa.faa")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            8
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/add_kraken_taxa_proteomes.log"
        shadow: "copy-minimal"
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
            12
        conda:
            os.path.join(ENVDIR, "IMP_protein_db.yaml")
        log:
            "logs/analysis/add_extra_proteomes.log"
        message:
            "Adding extra proteome(s) to the metaproteomics proteome ..."
        shadow: "copy-minimal"
        shell:
            """
            echo "Merging extra proteome(s) with metaproteomics proteome ..." >> {log}
            (
                cat {input} 
                for file in {params.extra_proteomes}; do
                    if [[ $file == *.gz ]]; then
                        zcat "$file"
                    else
                        cat "$file"
                    fi
                done
            ) > {output} 2>> {log}
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
        runtime = SHORTISH_RUNTIME,
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
    shadow: "copy-minimal"
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
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        sample = config['sample']
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_protein_db.yaml")
    log:
        "logs/analysis/rename_proteome_duplicates.log"
    shadow: "copy-minimal"
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
        runtime = SHORTISH_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        crap_fasta = config['proteomics']['crap_db']
    threads:
        getThreads(12)
    log:
        "logs/analysis/add_crap_to_fasta.log"
    shadow: "copy-minimal"
    shell:
        """
        echo "Adding CRAP proteins to the protein fasta ..." >> {log}
        cat {input.renamed_fasta} \
            {params.crap_fasta} 2>> {log} \
          | sed -e "s/ /__/g" 2>> {log} \
          | pigz -p {threads} -c > {output.final_fasta} 2>> {log}
        echo "Done." >> {log}
        """
