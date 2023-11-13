if config['proteomics']['insert_variants']:
    rule insertVariants:
        input:
            annotation="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
            assembly="Assembly/%s.assembly.merged.fa" %ASS,
            var=expand("Analysis/snps/{type}.variants.samtools.vcf.gz",type=TYPES)
        output:
            predictions="Analysis/annotation/proteomics.proteins.faa",
            peptideTab="Analysis/annotation/proteomics.pos.tsv",
            addInfo="Analysis/annotation/proteomics.info"
        resources: 
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            minPep = config['proteomics']['filter_N_peptides'],
            out = "Analysis/annotation/proteomics"
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        threads:
            1
        log:
            "logs/analysis_insertVariants.log"
        message:
            "insertVariants: inserting variants in protein predictions."
        shell:
            """
            export PATH=${{CONDA_PREFIX}}/bin:$PATH
            export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
            {SRCDIR}/variant_monsterProtein_anvioFilter_filter.pl -a {input.annotation} \
                                                                  -s {input.assembly} \
                                                                  -v {input.var} \
                                                                  -u {params.minPep} \
                                                                  -n {params.out} \
                                                                  -p \
                                                                  -anvio >> {log} 2>&1
            """
elif config['proteomics']['filter_N_peptides']:
    rule cutends:
        input:
            "Analysis/annotation/prokka.faa",
            "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
        output:
            "Analysis/annotation/proteomics.proteins.faa"
        resources: 
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            minPep = config['proteomics']['filter_N_peptides']
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        threads:
            1
        log:
            "logs/analysis_cutends.log"
        message:
            "cutends: removing non-tryptic peptides from protein predictions."
        shell:
            """
            export PATH=${{CONDA_PREFIX}}/bin:$PATH
            export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
            {SRCDIR}/trypsinStartEnd.filter.pl {input} \
                                               6 \
                                               {params.minPep} >> {output}
            """
else:
    rule renameProteomics:
        input:
            "Analysis/annotation/prokka.faa"
        output:
            "Analysis/annotation/proteomics.proteins.faa"
        resources: 
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_bioperl.yaml")
        threads:
            1
        log:
            "logs/analysis_renameProteomics.log"
        message:
            "renameProteomics: renaming protein predictions."
        shell:
            """
            export PATH=${{CONDA_PREFIX}}/bin:$PATH
            export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
            {SRCDIR}/rename4proteomics.pl {input} >> {output}
            """

if config['proteomics']['host_proteome'].split():
    rule mergeProteomicsHost:
        input:
            "Analysis/annotation/proteomics.proteins.faa"
        output:
            "Analysis/annotation/proteomics.final.faa"
        resources: 
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            host = config['proteomics']['host_proteome'].split()
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_fasta.yaml")
        log:
            "logs/analysis_mergeProteomicsHost.log"
        message:
            "mergeProteomicsHost: merging metaomics and host protein predictions."
        shell:
            """
            cat {input} {params.host} | fasta_formatter -w 60 >> {output}
            """
else:
    rule formatProteomics:
        input:
            "Analysis/annotation/proteomics.proteins.faa"
        output:
            "Analysis/annotation/proteomics.final.faa"
        resources: 
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_fasta.yaml")
        log:
            "logs/analysis_formatProteomics.log"
        message:
            "formatProteomics."
        shell:
            """
            cat {input} | fasta_formatter -w 60 >> {output}
            """

