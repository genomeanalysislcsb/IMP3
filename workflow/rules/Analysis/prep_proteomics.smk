if config["proteomics"]["insert_variants"]:
    rule analysis_insertVariants:
        input:
            annotation="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
            assembly="Assembly/%s.assembly.merged.fa" %ASS,
            var=expand("Analysis/snps/{type}.variants.samtools.vcf.gz",type=TYPES)
        output:
            faa="Analysis/annotation/proteomics.proteins.faa",
            tsv="Analysis/annotation/proteomics.pos.tsv",
            info="Analysis/annotation/proteomics.info",
        log:
            "logs/analysis_insertVariants.log"
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
            "insertVariants: inserting variants in proteins: {input}"
        shell:
            "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: avoid hardcoded paths
            "{params.script} -a {input.annotation} -s {input.assembly} -v {input.var} -u {params.minPep} -n $(dirname {output.faa}) -p -anvio &> {log}"
elif config["proteomics"]["filter_n_peptides"]:
    rule analysis_cutEnds:
        input:
            "Analysis/annotation/prokka.faa",
            "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
        output:
            "Analysis/annotation/proteomics.proteins.faa"
        log:
            "logs/analysis_cutends.log"
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
            "cutends: cut non-tryptic peptides from proteins: {input}"
        shell:
            "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: avoid hardcoded paths
            "{params.script} {input} 6 {params.minPep} > {output} 2> {log}" # TODO: why ">> {output}" ??? 
else:
    rule analysis_rename4proteomics:
        input:
            "Analysis/annotation/prokka.faa"
        output:
            "Analysis/annotation/proteomics.proteins.faa"
        log:
            "logs/analysis_renameProteomics.log"
        params:
            script=os.path.join(SRCDIR, "rename4proteomics.pl")
        threads: 1
        resources: 
            runtime="04:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        message:
            "renameProteomics: rename proteins: {input}"
        shell:
            "export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2 && " # TODO: avoid hardcoded paths
            "{params.script} {input} > {output} 2> {log}" # TODO: why ">> {output}" ??? 

if config["proteomics"]["host_proteome"]:
    rule analysis_mergeProteomicsHost:
        input:
            metap="Analysis/annotation/proteomics.proteins.faa",
            hostp=config["proteomics"]["host_proteome"], # TODO: moved from params to input -> would that be okay?
        output:
            "Analysis/annotation/proteomics.final.faa"
        log:
            "logs/analysis_mergeProteomicsHost.log"
        threads: 1
        resources: 
            runtime="04:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_fasta.yaml")
        message:
            "mergeProteomicsHost: merging metaP and host proteins: {input}"
        shell:
            "cat {input} | fasta_formatter -w 60 > {output} 2> {log}" # TODO: why ">> {output}" ??? 
else:
    rule formatProteomics:
        input:
            "Analysis/annotation/proteomics.proteins.faa"
        output:
            "Analysis/annotation/proteomics.final.faa"
        log:
            "logs/analysis_formatProteomics.log"
        threads: 1
        resources: 
            runtime="04:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_fasta.yaml")
        message:
            "formatProteomics: {input}"
        shell:
            "cat {input} | fasta_formatter -w 60 > {output} 2> {log}" # TODO: why ">> {output}" ??? 
