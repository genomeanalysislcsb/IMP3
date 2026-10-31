rule get_proteins_DAStool:
    input:
        'Assembly/%s.assembly.merged.fa' % ASS,
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Binning/prokka.renamed.faa"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        ENVDIR + "/IMP_bioperl.yaml"
    log:
        "logs/binning/translateProteins.log"
    message:
        "get_proteins_DAStool: Translating sequences based on GFF - code 11."
    shadow: "copy-minimal"
    shell:
        """
        export PATH=${{CONDA_PREFIX}}/bin:$PATH
        export PERL5LIB=${{CONDA_PREFIX}}/lib/perl5/site_perl
        {SRCDIR}/translateBacFromGFF.pl  {input[0]} {input[1]} >> {output} 2> {log}
        """

