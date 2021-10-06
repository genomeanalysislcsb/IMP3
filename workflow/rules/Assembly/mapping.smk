rule assembly_map_reads:
    input:
        r1="Preprocessing/{omic}.r1.preprocessed.fq",
        r2="Preprocessing/{omic}.r2.preprocessed.fq",
        se="Preprocessing/{omic}.se.preprocessed.fq",
        ref="Assembly/{asm}.assembly.merged.fa",
        idx=expand("Assembly/{{asm}}.assembly.merged.fa.{ext}", ext=get_bwa_idx_ext())
    output:
        temp("Assembly/{omic}.reads.bam")
    log:
        "logs/assembly.mapping.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
        asm=get_omics_assembly(),
    params:
        header=get_mapping_header,
    threads: get_cores_big()
    resources:
        runtime="12:00:00",
        mem=get_mem_big(get_cores_big()),
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Assembly: map reads on merged assembly: {input}"
    shell:
        """
        (samtools merge --threads {threads} -f {output} \
         <(bwa mem -v 1 -t {threads} -M -R \"{params.header}\" {input.ref} {input.r1} {input.r2} | samtools view --threads {threads} -bS -) \
         <(bwa mem -v 1 -t {threads} -M -R \"{params.header}\" {input.ref} {input.se} | samtools view --threads {threads} -bS -)) 2> {log}
        """

# TODO: should input/output be more generic ???
rule assembly_sort_bam:
    input:
        "Assembly/{omic}.reads.bam"
    output:
        "Assembly/{omic}.reads.sorted.bam"
    log:
        "logs/assembly.samtools.sort.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
    params:
        mem=get_mem_samtools_sort(),
    threads: get_cores_big()
    resources:
        runtime="12:00:00",
        mem=get_mem_big(get_cores_big()),
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Assembly: sort BAM: {input}"
    shell:
        "samtools sort --threads {threads} -m {params.mem} {input} > {output} 2> {log}"

rule get_unmapped_reads:
    input:
        r1="Preprocessing/{omic}.r1.preprocessed.fq",
        r2="Preprocessing/{omic}.r2.preprocessed.fq",
        se="Preprocessing/{omic}.se.preprocessed.fq",
        fa="Assembly/intermediary/{omic}.%s_preprocessed.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        idx=expand(
            "Assembly/intermediary/{{omic}}.%s_preprocessed.1.fa.{ext}" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            ext=get_bwa_idx_ext()
        ),
    output:
        r1="Assembly/intermediary/{omic}.r1.unmapped.fq",
        r2="Assembly/intermediary/{omic}.r2.unmapped.fq",
        se="Assembly/intermediary/{omic}.se.unmapped.fq",
        bam_tmp=temp("Assembly/intermediary/{omic}.pe.bam")
    log:
        "logs/assembly_extract_unmapped_from_%s.{omic}.log" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
    wildcard_constraints:
        omic="|".join(get_omics()),
        asm=get_omics_assembly(),
    params:
        mem=get_mem_samtools_sort(),
    threads: get_cores_big()
    resources:
        runtime="12:00:00",
        mem=get_mem_big(get_cores_big()),
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Extract reads from {input.r1}, {input.r2}, {input.se} not mapping to {input.fa}"
    shell: # TODO
        """
        bwa mem -v 1 -t {threads} {input.fa} {input.r1} {input.r2} 2>> {log} | \
        samtools view --threads {threads} -bS - > {output.bam_tmp} 2>> {log}
        samtools merge --threads {threads} -u - \
        <(samtools view --threads {threads} -u  -f 4 -F 264 {output.bam_tmp} 2>> {log}) \
        <(samtools view --threads {threads} -u -f 8 -F 260 {output.bam_tmp} 2>> {log}) \
        <(samtools view --threads {threads} -u -f 12 -F 256 {output.bam_tmp} 2>> {log}) 2>> {log} | \
        samtools view --threads {threads} -bF 0x800 - 2> {log} | samtools sort --threads {threads} -m {params.mem} -n - 2>> {log}| \
        bamToFastq -i stdin -fq {output.r1} -fq2 {output.r2} >> {log} 2>&1

        if [[ -s {input.se} ]]; then
            bwa mem -v 1 -t {threads} {input.fa} {input.se} 2>> {log} | \
            samtools view --threads {threads} -bS - 2>> {log} | samtools view --threads {threads} -uf 4 - 2>> {log} | \
            bamToFastq -i stdin -fq {output.se} >> {log} 2>&1
        else
            echo "There are no singletons reads. {input.se} is empty. Generating empty {output.se}" >> {log}
            touch {output.se}
        fi
        """

use rule get_unmapped_reads as get_unmapped_reads_hybrid with:
    input:
        r1="Preprocessing/{omic}.r1.preprocessed.fq",
        r2="Preprocessing/{omic}.r2.preprocessed.fq",
        se="Preprocessing/{omic}.se.preprocessed.fq",
        fa="Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        idx=expand(
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa.{ext}" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            ext=get_bwa_idx_ext()
        ),
    output:
        r1="Assembly/intermediary/mgmt.r1.%s_hybrid.{omic}.unmapped.fq" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        r2="Assembly/intermediary/mgmt.r2.%s_hybrid.{omic}.unmapped.fq" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        se="Assembly/intermediary/mgmt.se.%s_hybrid.{omic}.unmapped.fq" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        bam_tmp=temp("Assembly/intermediary/mgmt.pe.%s_hybrid.{omic}.bam" % get_assembler()), # TODO: rm assembler from the name or add as wildcard),
    log:
        "logs/assembly_extract_unmapped_from_hybrid_assembly.%s.{omic}.log" % get_assembler(), # TODO: rm assembler from the name or add as wildcard

#     use rule get_unmapped_reads as get_unmapped_reads_single_omic with:
#         input:
#             r1="Preprocessing/{omic}.r1.preprocessed.fq",
#             r2="Preprocessing/{omic}.r2.preprocessed.fq",
#             se="Preprocessing/{omic}.se.preprocessed.fq",
#             fa="Assembly/intermediary/{omic}.%s_preprocessed.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
#             idx=expand(
#                 "Assembly/intermediary/{{omic}}.%s_preprocessed.1.fa.{ext}" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
#                 ext=get_bwa_idx_ext()
#             ),
#         output:
#             r1="Assembly/intermediary/{omic}.r1.unmapped.fq",
#             r2="Assembly/intermediary/{omic}.r2.unmapped.fq",
#             se="Assembly/intermediary/{omic}.se.unmapped.fq"
#         log: "logs/assembly_extract_unmapped.{omic}.log"
