SAMTOOLS_MEM = str(round(float(MAIN_BIGMEMCORE[:-1]) * 0.75 - 0.5)) + config['cluster_specs']['mem_unit']

EXTRACT_UNMAPPED_SHELL = """
                         export TMPDIR="{params.tmp_dir}"
                         TMP_FILE=$(mktemp "alignment_XXXXXX.bam")

                         echo "Getting unmapped PE reads." >> {log} 2>&1
                         bwa-mem2 mem -v 1 \
                                      -t {threads} \
                                      {input[3]} \
                                      {input[0]} \
                                      {input[1]} 2>> {log} | \
                           samtools view --threads {threads} \
                                         --bam \
                                         - > ${{TMP_FILE}} 2>> {log}

                         samtools merge --threads {threads} \
                                        -u \
                                        - \
                                        <(samtools view --threads {threads} \
                                                        --uncompressed \
                                                        --require-flags 4 \
                                                        --exclude-flags 264 \
                                                        ${{TMP_FILE}} 2>> {log}) \
                                        <(samtools view --threads {threads} \
                                                        --uncompressed \
                                                        --require-flags 8 \
                                                        --exclude-flags 260 \
                                                        ${{TMP_FILE}} 2>> {log}) \
                                        <(samtools view --threads {threads} \
                                                        --uncompressed \
                                                        --require-flags 12 \
                                                        --exclude-flags 256 \
                                                        ${{TMP_FILE}} 2>> {log}) 2>> {log} | \
                           samtools view --threads {threads} \
                                         --bam \
                                         --exclude-flags 2048  \
                                         - 2>> {log} | \
                           samtools sort --threads {threads} \
                                         -m {SAMTOOLS_MEM} \
                                         -n \
                                         - 2>> {log} | \
                           bamToFastq -i stdin \
                                      -fq {output[0]} \
                                      -fq2 {output[1]} >> {log} 2>&1
                         rm -f ${{TMP_FILE}}

                         if [[ -s {input[2]} ]]; then
                           echo "Getting unmapped SE reads." >> {log} 2>&1

                           # atm piping directly from bwa --> samtools view --> bamToFastq fails for this part specifically ??
                           TMP_FILE_SE=$(mktemp "ex_unm_se_XXXXX.sam")
                           TMP_FILE_SE2=$(mktemp "ex_unm_se_XXXXX.bam")
                            
                           bwa-mem2 mem -v 1 \
                                        -t {threads} \
                                        {input[3]} \
                                        {input[2]} > ${{TMP_FILE_SE}} 2>> {log}
                            
                           samtools view --threads {threads} \
                                         --bam \
                                         --uncompressed \
                                         --require-flags 4 \
                                         ${{TMP_FILE_SE}} > ${{TMP_FILE_SE2}} 2>> {log}
                           
                           bamToFastq -i ${{TMP_FILE_SE2}} \
                                      -fq {output[2]} &>> {log}

                            rm -f ${{TMP_FILE_SE}} ${{TMP_FILE_SE2}}
                         else
                           echo "There are no singletons reads. {input[2]} is empty. Generating empty {output[2]}" >> {log}
                           touch {output[2]}
                         fi
                         """

rule extract_unmapped_from_megahit:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.0123',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.amb',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.ann',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.bwt.2bit.64',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.pac'
    output:
        'Assembly/intermediary/{type}.r1.unmapped.fq',
        'Assembly/intermediary/{type}.r2.unmapped.fq',
        'Assembly/intermediary/{type}.se.unmapped.fq'
    resources:
        runtime = "08:00:00",
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        tmp_dir=TMPDIR
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    log:
        "logs/assembly_extract_unmapped_from_megahit.{type}.log"
    message:
        "extract_unmapped_from_megahit: Extracting unmapped {wildcards.type} reads from megahit assembly."
    shell:
        EXTRACT_UNMAPPED_SHELL

rule extract_unmapped_from_hybrid_assembly:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.0123',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.amb',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.ann',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.bwt.2bit.64',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.pac'
    output:
        'Assembly/intermediary/mgmt.r1.{assembler}_hybrid.{type}.unmapped.fq',
        'Assembly/intermediary/mgmt.r2.{assembler}_hybrid.{type}.unmapped.fq',
        'Assembly/intermediary/mgmt.se.{assembler}_hybrid.{type}.unmapped.fq'
    resources:
        runtime = "12:00:00",
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        tmp_dir=TMPDIR
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    log:
        "logs/assembly_extract_unmapped_from_hybrid_assembly.{assembler}.{type}.log"
    message:
        "extract_unmapped_from_hybrid_assembly: Extracting unmapped {wildcards.type} reads from hybrid {wildcards.assembler} assembly."
    shell:
        EXTRACT_UNMAPPED_SHELL
