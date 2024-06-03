import math

SAMTOOLS_MEM = str(round(float(MAIN_BIGMEMCORE[:-1]) * 0.75 - 0.5)) + config['cluster_specs']['mem_unit']

rule mapping_on_assembly:
    input:
        r_1='Preprocessing/{type}.r1.preprocessed.fq',
        r_2='Preprocessing/{type}.r2.preprocessed.fq',
        r_se='Preprocessing/{type}.se.preprocessed.fq',
        assembly='Assembly/%s.assembly.merged.fa' % ASS,
        assembly_0123="Assembly/%s.assembly.merged.fa.0123" % ASS,
        assembly_amb="Assembly/%s.assembly.merged.fa.amb" % ASS,
        assembly_ann="Assembly/%s.assembly.merged.fa.ann" % ASS,
        assembly_bwt="Assembly/%s.assembly.merged.fa.bwt.2bit.64" % ASS,
        assembly_pac="Assembly/%s.assembly.merged.fa.pac" % ASS,
    output:
        map_pair = 'Assembly/{type}.reads.paired.sorted.bam',
        map_se = 'Assembly/{type}.reads.single.sorted.bam',
        map_merge = 'Assembly/{type}.reads.sorted.bam'
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
        "logs/assembly/mapping_on_assembly.{type}.log"
    message:
        "mapping_on_assembly: Mapping {wildcards.type} reads on merged assembly."
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
        echo "#####################################################################################" >> {log} 2>&1
        echo "Done." >> {log}
        """

rule index_merged_bam:
    input:
        bam = 'Assembly/{type}.reads.sorted.bam'
    output:
        bai = 'Assembly/{type}.reads.sorted.bam.bai'
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        # Ensure at least 50gb of mem
        max(int(math.ceil(50 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    log:
        "logs/assembly/index_bam.{type}.log"
    message:
        "index_bam: Indexing {wildcards.type} bam."
    shell:
        """
        echo "Indexing {input} ..." >> {log}
        samtools index {input} >> {log} 2>&1

        # Check output and print name
        if [ -f {output} ]; then
            echo "Indexing done. Index at: {output}." >> {log}
        else
            echo "Indexing failed." >> {log}
            exit 1
        fi

        echo "Done." >> {log}
        """

use rule index_merged_bam as rule_index_paired_bam with:
    input:
        bam = 'Assembly/{type}.reads.paired.sorted.bam'
    output:
        bai = 'Assembly/{type}.reads.paired.sorted.bam.bai'
    log:
        "logs/assembly/index_bam.{type}.paired.log"
    message:
        "index_bam: Indexing {wildcards.type} paired bam."

use rule index_merged_bam as rule_index_se_bam with:
    input:
        bam = 'Assembly/{type}.reads.single.sorted.bam'
    output:
        bai = 'Assembly/{type}.reads.single.sorted.bam.bai'
    log:
        "logs/assembly/index_bam.{type}.single.log"
    message:
        "index_bam: Indexing {wildcards.type} single bam."
