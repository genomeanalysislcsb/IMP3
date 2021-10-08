rule analysis_annotate:
    input:
        "Assembly/%s.assembly.merged.fa" % get_omics_assembly()
    output:
        fasta=expand("Analysis/annotation/prokka.{ext}", ext=get_prokka_ext("fasta")),
        gff=temp("Analysis/annotation/prokka.gff"),
    log:
        "logs/analysis_annotate.log"
    params:
        script=os.path.join(SRCDIR, "prokkaC")
    threads: 12
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(12),
    conda:
        os.path.join(ENVDIR + "/IMP_annotation.yaml")
    message:
        "Annotation: {input}"
    shell:
        "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: hardcoded paths should be avoided
        "export LC_ALL=en_US.utf-8 && "
        "if [ ! -f $CONDA_PREFIX/db/hmm/HAMAP.hmm.h3m ]; then {params.script} --dbdir $CONDA_PREFIX/db --setupdb; fi && " 
        "{params.script} --dbdir $CONDA_PREFIX/db --force --outdir $(dirname {output.gff}) --prefix prokka --noanno --cpus {threads} --metagenome {input} &> {log}"

rule analysis_filt_gff:
    input:
        "Analysis/annotation/prokka.gff"
    output:
        gff_gz="Analysis/annotation/prokka.gff.gz",
        gff="Analysis/annotation/annotation.filt.gff",
    threads: 1
    resources:
        runtime="01:00:00",
        mem=get_mem_norm(1),
    message:
        "Annotation: reformat GFF: {input}"
    shell:
        # remove prokka's header and FASTA part at the end of the file
        # TODO: simplify ???
        """
        LN=`grep -Hn "^>" {input} | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN {input} | grep -v "^#" | sort | uniq | grep -v "^==" > {output.gff}
        gzip -c {input} > {output.gff_gz}
        """

rule analysis_hmmsearch:
    input:
        "Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/prokka.faa.{db}.hmmscan"
    log:
        "logs/analysis_hmmer.{db}.log"
    wildcard_constraints:
        db="|".join(get_hmms())
    params: 
        params=get_hmm_cutoff,
        db_path=get_hmm_path,
    threads: 12
    resources:
        runtime="01:00:00",
        mem=get_mem_norm(12),
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "HMM search: {input}, DB={wildcards.db}"
    shell:
        "hmmsearch --cpu {threads} --noali --notextw {params.params} --tblout {output} {params.db_path}/*.hmm {input} > /dev/null 2> {log}"

rule analysis_hmmsearch_gff:
    input:
        hits="Analysis/annotation/intermediary/prokka.faa.{db}.hmmscan",
        gff="Analysis/annotation/annotation.filt.gff",
        faa="Analysis/annotation/prokka.faa"
    output:
        "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    log:
        "logs/analysis_makegff.{db}.log"
    wildcard_constraints:
        db="|".join(get_hmms())
    params:
        script=os.path.join(SRCDIR, "hmmscan_addBest2gff.pl"),
        # "{db}", # TODO: wildcards.db ???
        params=get_hmm_trim,
    threads: 1
    resources:
        runtime="4:00:00",
        mem=get_mem_norm(1),
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "HMM search GFF: {input}"
    shell:
        "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: hardcoded paths should be avoided
        "{params.script} -i {input.hits} -a {input.gff} -n {wildcards.db} {params.params} -o {output} -g $(grep \">\" {input.faa} | wc -l) &> {log}"

rule analysis_merge_gffs:
    input:
        gff="Analysis/annotation/annotation.filt.gff",
        hits=expand("Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff", db=get_hmms())
    output:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    log:
        "logs/analysis_mergegff.log"
    params:
        script=os.path.join(SRCDIR, "mergegffs.pl")
    threads: 1
    resources:
        runtime="8:00:00",
        mem=get_mem_norm(1),
    message:
        "Merge GFFs: {input}"
    shell:
        "{params.script} {output} {input} &> {log}"

rule analysis_call_variants:
    input:
        fa="Assembly/%s.assembly.merged.fa" % get_omics_assembly(),
        bam="Assembly/{omic}.reads.sorted.bam"
    output:
        vcf="Analysis/snps/{omic}.variants.samtools.vcf.gz",
        vcf_mpu=temp("Analysis/snps/{omic}.variants.samtools.mpu.vcf"),
        # vcf_plt=temp("Analysis/snps/{omic}.variants.samtools.plt.vcf"), # TODO: purpose of this tmp file (apparently it is not used) ???
    log:
        "logs/analysis_variant_calling.{omic}.log"
    threads: 1
    resources:
        runtime="48:00:00",
        mem=get_mem_norm(1)
    conda:
        os.path.join(ENVDIR, "IMP_mapping.yaml")
    message:
        "Call variants: Calling variants using bcftools for {wildcards.omic} reads."
    shell:
        "bcftools mpileup -A -Ou -f {input.fa} {input.bam} 2>> {log} | bcftools call -vm -Ov > {output.vcf_mpu} 2> {log} && "
        "bgzip -c {output.vcf_mpu} > {output.vcf} 2>> {log} && "
        "tabix -f -p vcf {output.vcf} &>> {log}"
