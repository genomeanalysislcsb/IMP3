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
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            getThreads(MAIN_BIGCORENO, bigmem=True)
        params:
            # sam_header = f'@RG\tID:{SAMPLE}\tSM:mg\tPL:platform\tLB:library'
            sample = SAMPLE
        conda:
            ENVDIR + "/IMP_mapping.yaml"
        log:
            "logs/assembly/mapping_on_assembly.mg.log"
        message:
            "mapping_on_assembly: Mapping mg reads on merged assembly."
        shell:
            """
            # map paired and sort
            echo "Mapping PE reads and sorting them." >> {log} 2>&1
            bwa-mem2 mem -v 1 \
                         -t {threads} \
                         -M \
                         -R '@RG\tID:{params.sample}\tSM:mg\tPL:platform\tLB:library' \
                         {input.assembly} \
                         {input.r_1} \
                         {input.r_2} 2>> {log} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_pair} 2>> {log}

            echo "#####################################################################################" >> {log} 2>&1

            # map se and sort
            echo "Mapping SE reads and sorting them." >> {log} 2>&1
            bwa-mem2 mem -v 1 \
                         -t {threads} \
                         -M \
                         -R '@RG\tID:{params.sample}\tSM:mg\tPL:platform\tLB:library' \
                         {input.assembly} \
                         {input.r_se} 2>> {log} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_se} 2>> {log}
            
            echo "#####################################################################################" >> {log} 2>&1

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
            temp("Assembly/Genome"),
            temp("Assembly/genomeParameters.txt"),
            temp("Assembly/chrStart.txt"),
            temp("Assembly/chrName.txt"),
            temp("Assembly/chrNameLength.txt"),
            temp("Assembly/chrLength.txt")
        resources:
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            # Ensure at least 75gb of mem
            max(int(math.ceil(75 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
        params:
            fasta_file = lambda wildcards: "_".join(f'Assembly/{wildcards.get("ASS")}.assembly.merged.fa'.split("/")[-2:]),
            max_mem_bytes = int(75e9)
        conda:
            os.path.join(ENVDIR, "IMP_mapping.yaml")
        log:
            "logs/assembly/index_assembly_STAR.log"
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
            unmap_pair_1 = 'Assembly/mt.reads.paired.sorted.unmapped_1.bam',
            unmap_pair_2 = 'Assembly/mt.reads.paired.sorted.unmapped_2.bam',
            map_se = 'Assembly/mt.reads.single.sorted.bam',
            unmap_se = 'Assembly/mt.reads.single.sorted.unmapped.bam',
            map_merge = 'Assembly/mt.reads.sorted.bam'
        resources:
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            getThreads(128, bigmem=True)
        params:
            # sam_header = f'@RG\tID:{SAMPLE}\tSM:mt'
            sample = SAMPLE,
            out_dir = 'Assembly',
            # Calc max amount of mem availbale (in bytes)
            max_mem = int(round(MAIN_BIGMEMCORE_MB * getThreads(128, bigmem=True)) * 1e6)
        conda:
            ENVDIR + "/IMP_mapping.yaml"
        log:
            "logs/assembly/mapping_on_assembly.mt.log"
        message:
            "mapping_on_assembly: Mapping mt reads on merged assembly."
        shell:
            """
            # Remove leftover tmp files
            rm -rf {params.out_dir}/_STARtmp {params.out_dir}/Log.out {params.out_dir}/Log.progress.out  \
                   {params.out_dir}/Log.final.out {params.out_dir}/SJ.out.tab > {log} 2>&1
            
            # Increase max number of open file descriptors for STAR
            ulimit -n 100000

            # Map paired
            echo "Mapping PE reads and sorting them ..." >> {log} 2>&1
            # --seedPerWindowNmax 30
            STAR --runThreadN {threads} \
                 --runMode alignReads \
                 --genomeDir Assembly \
                 --outFilterMultimapNmax 100 \
                 --winAnchorMultimapNmax 150 \
                 --alignIntronMax 1 \
                 --outFileNamePrefix {params.out_dir}/ \
                 --outSAMtype BAM SortedByCoordinate \
                 --limitBAMsortRAM {params.max_mem} \
                 --outSAMheaderHD '@RG\tID:{params.sample}\tSM:mt\tPL:platform\tLB:library' \
                 --outBAMcompression -1 \
                 --outReadsUnmapped Fastx \
                 --outFilterScoreMinOverLread 0.33 \
                 --outFilterMatchNminOverLread 0.33 \
                 --readFilesIn {input.r_1} {input.r_2} >> {log} 2>&1

            # Move files to output
            mv {params.out_dir}/Aligned.sortedByCoord.out.bam {output.map_pair} >> {log} 2>&1
            mv {params.out_dir}/Unmapped.out.mate1 {output.unmap_pair_1} >> {log} 2>&1
            mv {params.out_dir}/Unmapped.out.mate2 {output.unmap_pair_2} >> {log} 2>&1

            # Remove tmp files
            rm -rf {params.out_dir}/_STARtmp {params.out_dir}/Log.out {params.out_dir}/Log.progress.out  \
                   {params.out_dir}/Log.final.out {params.out_dir}/SJ.out.tab >> {log} 2>&1

            # Map single
            echo "Mapping SE reads and sorting them ..." >> {log} 2>&1
            # --seedPerWindowNmax 30
            STAR --runThreadN {threads} \
                 --runMode alignReads \
                 --genomeDir Assembly \
                 --outFilterMultimapNmax 100 \
                 --winAnchorMultimapNmax 150 \
                 --alignIntronMax 1 \
                 --outFileNamePrefix {params.out_dir}/ \
                 --outSAMtype BAM SortedByCoordinate \
                 --limitBAMsortRAM {params.max_mem} \
                 --outSAMheaderHD '@RG\tID:{params.sample}\tSM:mt\tPL:platform\tLB:library' \
                 --outBAMcompression -1 \
                 --outReadsUnmapped Fastx \
                 --outFilterScoreMinOverLread 0.33 \
                 --outFilterMatchNminOverLread 0.33 \
                 --readFilesIn {input.r_se} >> {log} 2>&1

            # Move files to output
            mv {params.out_dir}/Aligned.sortedByCoord.out.bam {output.map_se} >> {log} 2>&1
            mv {params.out_dir}/Unmapped.out.mate1 {output.unmap_se} >> {log} 2>&1

            # Remove tmp files
            rm -rf {params.out_dir}/_STARtmp {params.out_dir}/Log.out {params.out_dir}/Log.progress.out  \
                   {params.out_dir}/Log.final.out {params.out_dir}/SJ.out.tab >> {log} 2>&1

            # Merge paired and se and sort
            echo "Merging mapping PE and SE reads and sorting them ..." >> {log} 2>&1
            samtools merge --threads {threads} \
                           -o - \
                          {output.map_pair} \
                          {output.map_se} \
              | samtools sort --threads {threads} \
                              -m {SAMTOOLS_MEM} \
                              - > {output.map_merge} 2>> {log}
            echo "Done." >> {log}
            """

rule index_bam:
    input:
        '{sorted_bam}'
    output:
        '{sorted_bam}.bai'
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        # Ensure at least 50gb of mem
        max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    params:
        out_file_name = lambda wildcards: wildcards.get("sorted_bam").split("/")[-1].replace(".bam", "")
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    # log:
    #     f"logs/assembly/index.{params.out_file_name}.log"
    message:
        "index_bam: Indexing {wildcards.sorted_bam}."
    shell:
        """
        log="logs/assembly/index.{params.out_file_name}.log"

        # Ensure log dir exists
        mkdir -p logs/assembly

        samtools index {input} > ${{log}} 2>&1
        """
