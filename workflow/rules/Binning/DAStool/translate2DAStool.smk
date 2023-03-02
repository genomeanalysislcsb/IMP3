rule get_proteins_DAStool:
    input:
        'Assembly/%s.assembly.merged.fa' % ASS,
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Binning/prokka.renamed.faa"
    resources:
        runtime = "02:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    conda:
        ENVDIR + "/IMP_bioperl.yaml"
    log:
        "logs/binning_translateProteins.log"
    message:
        "get_proteins_DAStool: Translating sequences based on GFF - code 11."
    shell:
        """
        export PATH=${{CONDA_PREFIX}}/bin:$PATH
        export PERL5LIB=$CONDA_PREFIX/lib/perl5
        {SRCDIR}/translateBacFromGFF.pl  {input[0]} {input[1]} >> {output} 2> {log}
        """

