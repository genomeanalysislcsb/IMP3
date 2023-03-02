rule assembly_contig_length:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Stats/%s/%s.assembly.length.txt" % (ASS,ASS),
        "Stats/%s/%s.assembly.gc_content.txt" % (ASS,ASS)
    resources:
        runtime = "02:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis_assembly_contig_length.log"
    message:
        "assembly_contig_length: Getting data on assembly length and GC content."
    shell:
        """
        export PATH=${{CONDA_PREFIX}}/bin:$PATH
        export PERL5LIB=$CONDA_PREFIX/lib/perl5
        perl {SRCDIR}/fastaNamesSizes.pl {input} > {output[0]} 2>> {log}

        echo "Obtaining GC content"
        # AHB TMP_GC=$(mktemp --tmpdir={TMPDIR} -t "gc_out_XXXXXX.txt")
        TMP_GC=$(mktemp --tmpdir=Analysis -t "gc_out_XXXXXX.txt")
        perl {SRCDIR}/get_GC_content.pl {input} $TMP_GC >> {log} 2>&1

        # The program above provides a file gc_out.txt. This command cleans the output
        echo "Clean up output" >> {log}
        cut -f1,2 $TMP_GC | sed -e 's/>//g'> {output[1]}
        echo "Remove intermediate files" >> {log}
        rm $TMP_GC
        """
