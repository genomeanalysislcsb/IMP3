localrules: allmotu

rule allmotu:
    input:
        expand("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv", type=TYPES)
    output:
        touch("status/taxonomy_mOTUs.done")
    threads: 1
    message: "allmotu: mOTUs done?"

rule motus:
    input:
        "Preprocessing/{type}.r1.preprocessed.fq",
        "Preprocessing/{type}.r2.preprocessed.fq",
        "Preprocessing/{type}.se.preprocessed.fq"
    output:
        report("Analysis/taxonomy/mOTUs/{type}.mOTU.counts.tsv",category="Taxonomy")
    threads: getThreads(6)
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_taxonomy.yaml"
    log: "logs/analysis_motus.{type}.log"
    message: "motus: Running mOTUs on {wildcards.type} reads."
    shell:
        """
        a=`stat -c %A $CONDA_PREFIX/share/motus-*/db_mOTU/db_mOTU_versions`
        if [ $a == "-rw-------" ]
	then
	    chmod 666 $CONDA_PREFIX/share/motus-*/db_mOTU/*
	fi
        if [ -s {input[2]} ] && [ -s {input[0]} ] && [ -s {input[1]} ]
        then
            motus profile -f {input[0]} -r {input[1]} -s {input[2]} -A -c -t {threads} -o {output} >& {log}
        elif [ -s {input[2]} ]
        then
            motus profile -s {input[2]} -A -c -t {threads} -o {output} >& {log}
        else
            motus profile -f {input[0]} -r {input[1]} -A -c -t {threads} -o {output} >& {log}
        fi 
        """

