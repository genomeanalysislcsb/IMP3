import math

SAMTOOLS_MEM = str(round(float(MAIN_BIGMEMCORE[:-1]) * 0.75 - 0.5)) + config['cluster_specs']['mem_unit']

if 'mg' in TYPES:
    rule mg_mapping_on_assembly:
        input:
            r_1='Preprocessing/mg.r1.preprocessed.fq',
            r_2='Preprocessing/mg.r2.preprocessed.fq',
            r_se='Preprocessing/mg.se.preprocessed.fq',
            assembly='Assembly/%s.assembly.merged.fa' % ASS,
            assembly_0123="Assembly/%s.assembly.merged.fa.0123" % ASS,
            assembly_amb="Assembly/%s.assembly.merged.fa.amb" % ASS,
            assembly_ann="Assembly/%s.assembly.merged.fa.ann" % ASS,
            assembly_bwt="Assembly/%s.assembly.merged.fa.bwt.2bit.64" % ASS,
            assembly_pac="Assembly/%s.assembly.merged.fa.pac" % ASS,
        output:
            map_pair = 'Assembly/mg.reads.paired.sorted.bam',
            map_se = 'Assembly/mg.reads.single.sorted.bam',
            map_merge = 'Assembly/mg.reads.sorted.bam'
        resources:
            runtime = "12:00:00",
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            getThreads(MAIN_BIGCORENO, bigmem=True)
        params:
            sam_header = f'@RG\\tID:{SAMPLE}\\tSM:mg'
        conda:
            ENVDIR + "/IMP_mapping.yaml"
        log:
            "logs/assembly_mapping_on_assembly.mg.log"
        message:
            "mapping_on_assembly: Mapping mg reads on merged assembly."
        shell:
            """
            # map paired and sort
            echo "Mapping PE reads and sorting them." >> {log} 2>&1
            bwa-mem2 mem -v 1 \
                         -t {threads} \
                         -M \
                         -R {params.sam_header} \
                         {input.assembly} \
                         {input.r_1} \
                         {input.r_2} 2>> {log} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_pair} 2>> {log}

            echo "#####################################################################################" >> {log}

            # map se and sort
            echo "Mapping SE reads and sorting them." >> {log} 2>&
            bwa-mem2 mem -v 1 \
                         -t {threads} \
                         -M \
                         -R {params.sam_header} \
                         {input.assembly} \
                         {input.r_se} 2>> {log} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_se} 2>> {log}
            
            echo "#####################################################################################" >> {log}

            # merge paired and se and sort
            echo "Merging mapping PE and SE reads and sorting them." >> {log} 2>&1
            samtools merge --threads {threads} \
                           -o - \
                          {output.map_pair} \
                          {output.map_se} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_merge} 2>> {log}
            """

if 'mt' in TYPES:
    rule star_index:
        input:
            r_1='Preprocessing/mt.r1.preprocessed.fq',
            assembly='Assembly/%s.assembly.merged.fa' % ASS
        output:
            directory("Assembly/Genome"),
            "Assembly/genomeParameters.txt",
            "Assembly/chrStart.txt",
            "Assembly/chrName.txt",
            "Assembly/chrNameLength.txt",
            "Assembly/chrLength.txt"
        resources:
            runtime = "24:00:00",
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            # Ensure at least 20gb of mem
            max(int(math.ceil(20 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
        params:
            fasta_file = lambda wc: "_".join(f'Assembly/{wc.get("ASS")}.assembly.merged.fa'.split("/")[-2:]),
            max_mem_bytes = int(20e9)
        conda:
            os.path.join(ENVDIR, "IMP_mapping.yaml")
        log:
            "logs/assembly_index_assembly_STAR.log"
        message:
            "STAR_index: Indexing {params.fasta_file} for STAR."
        shell:
            """
            # Get stats
            r_max_l=$( seqkit stats -T {input.r_1} | cut -f 8 | tail -n +2 )
            assembly_stats=$( seqkit stats -T {input.assembly} | tail -n +2 )
            a_contigs=$( echo ${{assembly_stats}} | cut -d ' ' -f 4 )
            a_length=$( echo ${{assembly_stats}} | cut -d ' ' -f 5 )

            genomeChrBinNbits=$( python -c "import math; print( round(min(18, math.log2(max(${{a_length}}/${{a_contigs}}, ${{r_max_l}})))) )" )
            genomeSAindexNbases=$( python -c "import math; print( round(min(14, math.log2(${{a_length}}) / 2 - 1)))" )

            echo "r_max_l: ${{r_max_l}}" >> {log}
            echo "assembly_stats: ${{assembly_stats}}" >> {log}
            echo "genomeChrBinNbits: ${{genomeChrBinNbits}}" >> {log}
            echo "genomeSAindexNbases: ${{genomeSAindexNbases}}" >> {log}

            STAR --runThreadN {threads} \
                 --runMode genomeGenerate \
                 --genomeDir Assembly \
                 --genomeChrBinNbits ${{genomeChrBinNbits}} \
                 --genomeSAindexNbases 10 \
                 --limitGenomeGenerateRAM {params.max_mem_bytes} \
                 --genomeFastaFiles {input.assembly} >> {log} 2>&1
            """
    
    rule mt_mapping_on_assembly:
        input:
            r_1='Preprocessing/mt.r1.preprocessed.fq',
            r_2='Preprocessing/mt.r2.preprocessed.fq',
            r_se='Preprocessing/mt.se.preprocessed.fq',
            star_in_dir="Assembly/Genome",
            star_in_param="Assembly/genomeParameters.txt",
            star_in_chr_st="Assembly/chrStart.txt",
            star_in_chr_n="Assembly/chrName.txt",
            star_in_chr_nl="Assembly/chrNameLength.txt",
            star_in_chr_l="Assembly/chrLength.txt",
            assembly='Assembly/%s.assembly.merged.fa' % ASS
        output:
            map_pair = 'Assembly/mt.reads.paired.sorted.bam',
            map_se = 'Assembly/mt.reads.single.sorted.bam',
            map_merge = 'Assembly/mt.reads.sorted.bam'
        resources:
            runtime = "12:00:00",
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            getThreads(MAIN_BIGCORENO, bigmem=True)
        params:
            sam_header = f'@RG\\tID:{SAMPLE}\\tSM:mt'
        conda:
            ENVDIR + "/IMP_mapping.yaml"
        log:
            "logs/assembly_mapping_on_assembly.mt.log"
        message:
            "mapping_on_assembly: Mapping mt reads on merged assembly."
        shell:
            """
            # Map paired
            STAR --runThreadN {threads} \
                 --runMode alignReads \
                 --genomeDir Assembly \
                 --outFilterMultimapNmax 100 \
                 --winAnchorMultimapNmax 150 \
                 --alignIntronMax 1 \
                 --outFileNamePrefix {output.map_pair} \
                 --outSAMtype BAM SortedByCoordinate \
                 --outSAMheaderHD {params.sam_header} \
                 --outBAMcompression -1 \
                 --outReadsUnmapped Fastx \
                 --outFilterScoreMinOverLread 0.33 \
                 --outFilterMatchNminOverLread 0.33 \
                 --readFilesIn {input.r_1} {input.r_2} > {log} 2>&1

            # Map single
            STAR --runThreadN {threads} \
                 --runMode alignReads \
                 --genomeDir ${out_dir} \
                 --outFilterMultimapNmax 100 \
                 --winAnchorMultimapNmax 150 \
                 --alignIntronMax 1 \
                 --outFileNamePrefix {output.map_se} \ \
                 --outSAMtype BAM SortedByCoordinate \
                 --outSAMheaderHD ${{SAMHEADER}} \
                 --outBAMcompression -1 \
                 --outReadsUnmapped Fastx \
                 --outFilterScoreMinOverLread 0.33 \
                 --outFilterMatchNminOverLread 0.33 \
                 --readFilesIn {input.r_se} >> {log} 2>&1

            # merge paired and se and sort
            samtools merge --threads {threads} \
                           -o - \
                          {output.map_pair} \
                          {output.map_se} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_merge} 2>> {log} 
            """

rule index_bam:
    input:
        '{sorted_bam}'
    output:
        '{sorted_bam}.bai'
    resources:
        runtime = "04:00:00",
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        1
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    log:
        "logs/assembly_index.{sorted_bam}.log"
    message:
        "index_bam: Indexing {wildcards.sorted_bam}."
    shell:
        """
        samtools index {input} > {log} 2>&1
        """
