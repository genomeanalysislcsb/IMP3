rule maxbin_contig2bin:
    input:
        "Dereplication/MaxBin/{sample}/maxbin_res.log"
    output:
        "Dereplication/MaxBin/{sample}/maxbin_contig2bin.txt",
        "Dereplication/MaxBin/{sample}/scaffold2bin.tsv"
    threads: 1
    resources:
        runtime = "3:00:00",
        mem = MEMCORE
    log: "logs/dereplication_maxbin_contig2bin.{sample}.log"
    message: "maxbin_contig2bin: Getting contig to bin mapping."
    shell:
        """
        for f in $(ls Dereplication/MaxBin/{wildcards.sample}/*.fasta | cut -f4 -d "/"); do 
         sed -e "s/$/\t$f/g" <(grep "^>" Dereplication/MaxBin/{wildcards.sample}/$f | sed -e 's/>//g') >> {output[0]}
        done
        ln -s ../../../{output[0]} {output[1]}
        echo "Contig to bin mapping done" > {log}
        """

rule maxbin_get_depth:
    input:
        "Dereplication/cross_mapping/{sample}.contig_depth.{cross_sample}.txt",
        lambda wildcards: expand("{raw_directory}/Assembly/{ass}.assembly.merged.fa",
                raw_directory=samples.path[wildcards.sample],
                ass=samples.assembly_type[wildcards.sample])
    output:
        "Dereplication/cross_mapping/{sample}.contig_depth0.{cross_sample}.txt"
    threads: 1
    resources:
        runtime = "3:00:00",
        mem = MEMCORE
    message: "maxbin_get_depth: Getting depth of coverage by {wildcards.cross_sample} on {wildcards.sample} contigs in Maxbin format."
    shell:
        """
        cat {input[0]} \
         <(awk 'BEGIN {{OFS="\t"}}; {{print $0, 0}}' \
         <(diff --new-line-format="" --unchanged-line-format="" \
         <(grep "^>" {input[1]} |sed -e 's/>//g' | sort) \
         <(cut -f1 {input[0]} | sort))) > {output[0]}
        """


localrules: maxbin_prep

rule maxbin_prep:
    input:
        lambda wildcards: expand("Dereplication/cross_mapping/{{sample}}.contig_depth0.{cross_sample}.txt",
               cross_sample=samples.loc[[a and b for a,b in zip(samples.mg == 1,samples["sample"] != wildcards.sample)],"sample"]),
        lambda wildcards: expand("{raw_directory}/Stats/mg/{ass}.assembly.contig_depth.0.txt",
               raw_directory=samples.path[wildcards.sample],
               ass=samples.assembly_type[wildcards.sample])
    output:
        "Dereplication/cross_mapping/{sample}.contig_depth.0.list"
    message: "maxbin_prep: Collecting abundance for {wildcards.sample}."
    shell:
        """
        ls {input} >> {output}
        """
        
rule maxbin:
    input:
       lambda wildcards: expand("{raw_directory}/Assembly/{ass}.assembly.merged.fa",
                raw_directory=samples.path[wildcards.sample],
                ass=samples.assembly_type[wildcards.sample]),
       "Dereplication/cross_mapping/{sample}.contig_depth.0.list" 
    output:
       "Dereplication/MaxBin/{sample}/maxbin_res.log",
       "Dereplication/MaxBin/{sample}/maxbin_res.summary"
    threads: getThreads(2)
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/dereplication_maxbin.{sample}.log"
    message: "maxbin: Running MaxBin on {wildcards.sample}."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
	## Run MaxBin
        run_MaxBin.pl -contig {input[0]} \
         -abund_list {input[1]} \
         -out Dereplication/MaxBin/{wildcards.sample}/maxbin_res \
         -thread {threads} \
         -min_contig_length {config[dereplication][cross_mapping_rebinning][MaxBin][cutoff]} > {log} 2>&1
        """

