SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

EXTRACT_UNMAPPED_SHELL = """
 TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
 bwa mem -v 1 -t {threads} {input[3]} {input[0]} {input[1]} 2>> {log} | \
  samtools view --threads {threads} -bS - > $TMP_FILE 2>> {log}
 samtools merge --threads {threads} -u - \
  <(samtools view --threads {threads} -u  -f 4 -F 264 $TMP_FILE 2>> {log}) \
  <(samtools view --threads {threads} -u -f 8 -F 260 $TMP_FILE 2>> {log}) \
  <(samtools view --threads {threads} -u -f 12 -F 256 $TMP_FILE 2>> {log}) 2>> {log} | \
  samtools view --threads {threads} -bF 0x800 - 2> {log} | samtools sort --threads {threads} -m {SAMTOOLS_MEM} -n - 2>> {log}| \
  bamToFastq -i stdin -fq {output[0]} -fq2 {output[1]} >> {log} 2>&1

 if [[ -s {input[2]} ]]
 then
    bwa mem -v 1 -t {threads} {input[3]} {input[2]} 2>> {log} | \
     samtools view --threads {threads} -bS - 2>> {log} | samtools view --threads {threads} -uf 4 - 2>> {log} | \
     bamToFastq -i stdin -fq {output[2]} >> {log} 2>&1
 else
    echo "There are no singletons reads. {input[2]} is empty. Generating empty {output[2]}" >> {log}
    touch {output[2]}
 fi
 rm -rf  $TMP_FILE
"""

rule extract_unmapped_from_megahit:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.amb',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.bwt',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.pac',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.sa',
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa.ann',
    output:
        'Assembly/intermediary/{type}.r1.unmapped.fq',
        'Assembly/intermediary/{type}.r2.unmapped.fq',
        'Assembly/intermediary/{type}.se.unmapped.fq'
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: getThreads(BIGCORENO)
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_extract_unmapped_from_megahit.{type}.log"
    message: "extract_unmapped_from_megahit: Extracting unmapped {wildcards.type} reads from megahit assembly."
    shell:
        EXTRACT_UNMAPPED_SHELL

rule extract_unmapped_from_hybrid_assembly:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.amb',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.bwt',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.pac',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.sa',
        'Assembly/intermediary/mgmt.{assembler}_hybrid.1.fa.ann',
    output:
        'Assembly/intermediary/mgmt.r1.{assembler}_hybrid.{type}.unmapped.fq',
        'Assembly/intermediary/mgmt.r2.{assembler}_hybrid.{type}.unmapped.fq',
        'Assembly/intermediary/mgmt.se.{assembler}_hybrid.{type}.unmapped.fq'
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: getThreads(BIGCORENO)
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_extract_unmapped_from_hybrid_assembly.{assembler}.{type}.log"
    message: "extract_unmapped_from_hybrid_assembly: Extracting unmapped {wildcards.type} reads from hybrid {wildcards.assembler} assembly."
    shell:
        EXTRACT_UNMAPPED_SHELL
