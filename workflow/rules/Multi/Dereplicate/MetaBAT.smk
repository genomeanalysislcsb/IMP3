rule metabat:
    input:
        lambda wildcards: expand("{raw_directory}/Assembly/{ass}.assembly.merged.fa",
                raw_directory=samples.path[wildcards.sample],
                ass=samples.assembly_type[wildcards.sample]),
        "Dereplication/cross_mapping/{sample}.contig_depth.txt" 
    output:
        "Dereplication/MetaBAT/{sample}/metabat_res",
        "Dereplication/MetaBAT/{sample}/scaffold2bin.tsv"
    threads: getThreads(2)
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/Dereplication_MetaBAT.{sample}.log"
    message: "Running MetaBAT for {wildcards.sample}."
    shell:
        """
        ## Run MetaBat
        metabat2 -i {input[0]} \
          -a {input[1]} \
	     --saveCls \
         -o {output[0]} \
         -t {threads} \
         -m {config[dereplication][cross_mapping_rebinning][MetaBAT][cutoff]} > {log} 2>&1
	    ln -s ../../../{output[0]} {output[1]}
        """

rule metabat_prep:
    input:
        lambda wildcards: expand("Dereplication/cross_mapping/{cross_sample}.on.{{sample}}.sorted.bam",
                cross_sample = samples.loc[[a and b for a,b in zip(samples.mg == 1,samples["sample"] != wildcards.sample)],"sample"]),
        lambda wildcards: expand("{raw_directory}/Assembly/mg.reads.sorted.bam",
                raw_directory=samples.path[wildcards.sample])
    output:
        "Dereplication/cross_mapping/{sample}.contig_depth.txt"
    threads: 1
    resources:
        runtime = "24:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/binning_MetaBAT.{sample}.log"
    message: "Running MetaBAT."
    shell:
        """
	    jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """
