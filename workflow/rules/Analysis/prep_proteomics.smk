if config["proteomics"]["insert_variants"]:
    rule analysis_prepprot_insert_variants:
        input:
            annotation="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
            assembly="Assembly/%s.assembly.merged.fa" %ASS,
            var=expand("Analysis/snps/{type}.variants.samtools.vcf.gz",type=TYPES)
        output:
            faa="Analysis/annotation/proteomics.proteins.faa",
            tsv="Analysis/annotation/proteomics.pos.tsv",
            info="Analysis/annotation/proteomics.info",
        log:
            "logs/analysis_prepprot_insert_variants.log"
        benchmark:
            "benchmark/analysis_prepprot_insert_variants.tsv"
        params:
            script=os.path.join(SRCDIR, "variant_monsterProtein_anvioFilter_filter.pl"),
            minPep=config["proteomics"]["filter_n_peptides"],
        threads: 1
        resources: 
            runtime="48:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        message:
            "Prepare proteomics: insert variants in proteins: {input}"
        shell:
            "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: avoid hardcoded paths
            "{params.script} -a {input.annotation} -s {input.assembly} -v {input.var} -u {params.minPep} -n $(dirname {output.faa}) -p -anvio &> {log}"
elif config["proteomics"]["filter_n_peptides"]:
    rule analysis_prepprot_cut_ends:
        input:
            "Analysis/annotation/prokka.faa",
            "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
        output:
            "Analysis/annotation/proteomics.proteins.faa"
        log:
            "logs/analysis_prepprot_cut_ends.log"
        benchmark:
            "benchmark/analysis_prepprot_cut_ends.tsv"
        params:
            script=os.path.join(SRCDIR, "trypsinStartEnd.filter.pl"),
            minPep=config["proteomics"]["filter_n_peptides"],
        threads: 1
        resources: 
            runtime="04:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        message:
            "Prepare proteomics: cut non-tryptic peptides from proteins: {input}"
        shell:
            "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: avoid hardcoded paths
            "{params.script} {input} 6 {params.minPep} > {output} 2> {log}" # TODO: why ">> {output}" ??? 
else:
    rule analysis_prepprot_rename:
        input:
            "Analysis/annotation/prokka.faa"
        output:
            "Analysis/annotation/proteomics.proteins.faa"
        log:
            "logs/analysis_prepprot_rename.log"
        params:
            script=os.path.join(SRCDIR, "rename4proteomics.pl")
        threads: 1
        resources: 
            runtime="04:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        message:
            "Prepare proteomics: rename proteins: {input}"
        shell:
            "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: avoid hardcoded paths
            "{params.script} {input} > {output} 2> {log}" # TODO: why ">> {output}" ??? 

if config["proteomics"]["host_proteome"]:
    rule analysis_prepprot_add_host:
        input:
            metap="Analysis/annotation/proteomics.proteins.faa",
            hostp=config["proteomics"]["host_proteome"], # TODO: moved from params to input -> would that be okay?
        output:
            temp("Analysis/annotation/proteomics.proteins.host.faa")
        threads: 1
        resources: 
            runtime="01:00:00",
            mem=get_mem_norm,
        message:
            "Prepare proteomics: add host proteome: {input}"
        shell:
            "cat {input} > {output}" # TODO: why ">> {output}" ??? 

rule analysis_prepprot_reformat_faa:
    input:
        "Analysis/annotation/proteomics.proteins.host.faa" if has_host_proteome() else "Analysis/annotation/proteomics.proteins.faa"
    output:
        "Analysis/annotation/proteomics.final.faa"
    log:
        "logs/analysis_prepprot_reformat_faa.log"
    threads: 1
    resources: 
        runtime="04:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_fasta.yaml")
    message:
        "Prepare proteomics: reformat final FAA: {input}"
    shell:
        "cat {input} | fasta_formatter -w 60 > {output} 2> {log}" # TODO: why ">> {output}" ??? 
