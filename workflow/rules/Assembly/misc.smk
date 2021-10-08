rule bwa_idx:
    input:
        "{fasta}"
    output:
        expand("{{fasta}}.{ext}", ext=get_bwa_idx_ext())
    threads: 1
    resources:
        runtime="24:00:00",
        mem=get_mem_norm(1),
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "bwa index: {input}"
    shell:
        "bwa index {input}"

rule bam_idx:
    input:
        "{sorted_bam}"
    output:
        "{sorted_bam}.bai"
    log:
        "logs/samtools.idx_bam.{sorted_bam}.log"
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_big(1),
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "samtools: index BAM: {input}"
    shell:
        "samtools index {input} &> {log}"

rule assembly_idx:
    input:
        "Assembly/{asm}.assembly.merged.fa"
    output:
        "Assembly/{asm}.assembly.merged.fa.fai"
    log:
        "logs/assembly_idx_{asm}.log"
    wildcard_constraints:
        asm=get_omics_assembly(),
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1),
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Assembly index: {input}"
    shell:
        "samtools faidx {input} &> {log}"

rule assembly_fasta2bed:
    input:
        "Assembly/{asm}.assembly.merged.fa.fai"
    output:
        "Assembly/{asm}.assembly.merged.fa.bed3"
    wildcard_constraints:
        asm=get_omics_assembly(),
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1),
    message:
        "Assembly: FASTA to BED: {input}"
    shell:
        "cat {input} | awk '{{print $1 \"\t0\t\" $2}}' > {output}"

rule assembly_cov:
    input:
        bam="Assembly/{omic}.reads.sorted.bam",
        bed="Assembly/{asm}.assembly.merged.fa.bed3",
    output:
        "Stats/{omic}/{asm}.assembly.contig_coverage.txt"
    log:
        "logs/assembly_cov.{asm}.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
        asm=get_omics_assembly(),
    threads: 1
    resources:
        runtime="12:00:00",
        mem=get_mem_big(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Assembly: coverage: {input}"
    shell:
        "coverageBed -b {input.bam} -a {input.bed} -sorted > {output} 2> {log}"

rule assembly_depth:
    input:
        bam="Assembly/{omic}.reads.sorted.bam",
        fasta="Assembly/{asm}.assembly.merged.fa",
    output:
        depth="Stats/{omic}/{asm}.assembly.contig_depth.txt",
        depth_tmp=temp("Stats/{omic}/{asm}.assembly.contig_depth.txt.tmp"),
    log:
        "logs/assembly_mapping_stats.{asm}.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
        asm=get_omics_assembly(),
    params:
        script=os.path.join(SRCDIR, "calcAvgCoverage.pl"),
    threads: 1
    resources:
        runtime="12:00:00",
        mem=get_mem_big(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Assembly: contig depth: {input}"
    shell:
        "(genomeCoverageBed -ibam {input.bam} | grep -v \"genome\" > {output.depth_tmp} && "
        # method adapted and modified from CONCOCT
        "{params.script} {output.depth_tmp} {input.fasta} > {output.depth}) 2> {log}"

rule assembly_mapping_stats:
    input:
        "Assembly/{omic}.reads.sorted.bam"
    output:
        report("Stats/{omic}/{asm}.assembly.contig_flagstat.txt", category="Assembly")
    log:
        "logs/assembly_flagstat.{asm}.{omic}.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
        asm=get_omics_assembly(),
    threads: 1
    resources:
        runtime="12:00:00",
        mem=get_mem_big(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Assembly: flagstat: {input}"
    shell:
        "(samtools flagstat {input} | cut -f1 -d ' ' > {output}) 2> {log}"

# TODO: https://git-r3lab.uni.lu/IMP/imp3/-/issues/52
rule assembly_len:
    input:
        "Assembly/{asm}.assembly.merged.fa"
    output:
        "Stats/{asm}/{asm}.assembly.length.txt"
    log:
        "logs/assembly_len.{asm}.log"
    wildcard_constraints:
        asm=get_omics_assembly(),
    params:
        script=os.path.join(SRCDIR, "fastaNamesSizes.pl")
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "Assembly: contig length: {input}"
    shell:
        "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && perl {params.script} {input} > {output} 2> {log}"

# TODO: https://git-r3lab.uni.lu/IMP/imp3/-/issues/52
rule assembly_gc:
    input:
        "Assembly/{asm}.assembly.merged.fa"
    output:
        gc="Stats/{asm}/{asm}.assembly.gc_content.txt",
        tmp=temp("Stats/{asm}/{asm}.assembly.gc_content.txt.tmp"),
    log:
        "logs/analysis_assembly_contig_length.{asm}.log"
    wildcard_constraints:
        asm=get_omics_assembly(),
    params:
        script=os.path.join(SRCDIR, "get_GC_content.pl"),
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "Assembly: contig GC content: {input}"
    shell:
        "perl {params.script} {input} {output.tmp} &> {log} && cut -f1,2 {output.tmp} | sed -e 's/>//g' > {output.gc}"

# archive intermedite assembly files
rule archive_assembly:
    input:
        "status/assembly.ready" # see rule in modules/Assembly.smk
    output:
        "Assembly/intermediary.tar.gz"
    log:
        "logs/assembly_archive.log"
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1),
    params:
        folder="Assembly/intermediary/"
    message:
        "tar_intermediary: Compressing intermediary steps of assembly."
    shell:
       "tar czf {output} {params.folder} && rm -r {params.folder} &> {log}"
