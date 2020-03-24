import pandas as pd

localrules: grid_all

rule grid_all:
    input: 
        gather_bins
    output:
        "status/binning_GRiD.done"
    shell:
        """
        touch {output}
        """

if not "taxonomy" in IMP_STEPS:
    rule GRID_summarize:
        input:
            gather_bins
        output:
            report("Stats/%s/bins.tsv" %ASS,"Binning")
        params:
            runtime = "12:00:00",
            mem = MEMCORE
        log: "logs/binning_GRID_summarize.log"
        message: "GRID_summarize: Summarizing binning results."
        threads: 1
        conda: ENVDIR + "/IMP_binning.yaml"
        script:
            SRCDIR + "/summarize_GRID.R"
    
    rule bins_tar:
        input:
            dir="Binning/selected_DASTool_bins",
            ctrl="Stats/%s/bins.tsv" %ASS,
            bins=gather_bins
        output:
            "Binning/per_bin_results.tar.gz"
        threads: 1
        params:
            runtime = "8:00:00",
            mem = MEMCORE
        log: "logs/binning_bins_tar.log"
        message: "bins_tar: Compressing per-bin results of GRiD."
        shell:
            """
            tar cvzf {output} {input.bins} && rm -r {input.bins} >> {log} 2>&1
            """



checkpoint getBins:
    input: 
        "Binning/selected_DASTool_summary.txt"
    output:
        "status/binning_GRiD.prepped"
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
#    log: "logs/binning_prepareGRiD.log"
    message: "getBins: Copying DAStool bins to single directories."
    run:
        selected_bins=pd.read_table("Binning/selected_DASTool_summary.txt")['bin']
        for BIN in selected_bins:
            shell("mkdir -p Binning/selected_DASTool_bins/{BIN} && cp Binning/selected_DASTool_bins/{BIN}.contigs.fa Binning/selected_DASTool_bins/{BIN}/contigs.fa")
        shell("touch {output}")
    
rule grid:
    input:
       "Binning/selected_DASTool_bins/{clusterID}/contigs.fa",
       "Binning/selected_DASTool_bins/{clusterID}/reads.fastq"
    output:
        directory("Binning/selected_DASTool_bins/{clusterID}/grid")
    params:
        clusterInput="Binning/selected_DASTool_bins/",
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_grid.yaml"
    log: "logs/binning_grid.{clusterID}.log"
    message: "grid: Running GRiD on {wildcards.clusterID}."
    shell:
        """
        grid single -r {params.clusterInput}{wildcards.clusterID} -o {output} -g {input[0]} > {log} 2>&1
        """

rule prepare_gridBed:
    input:
        "Binning/selected_DASTool_bins/{clusterID}/contigs.fa"
    output:
        "Binning/selected_DASTool_bins/{clusterID}/contigs.bed"
    params:
        runtime = "2:00:00",
        mem = MEMCORE
    threads: 1
    log: "logs/binning_prepare_gridBed.{clusterID}.log"
    message: "prepare_gridBed: Preparing .bed file for GRiD - {wildcards.clusterID}."
    shell:
        """
        {SRCDIR}/fastaBed.pl {input} >> {output} 2> {log}
        """

rule prepare_gridReads:
    input:
        "Assembly/mg.reads.sorted.bam",
        "Binning/selected_DASTool_bins/{clusterID}/contigs.bed"
    output:
        "Binning/selected_DASTool_bins/{clusterID}/reads.fastq"
    params:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: 1
    log: "logs/binning_prepare_gridReads.{clusterID}.log"
    message: "prepare_gridReads: Preparing fastq file for GRiD - {wildcards.clusterID}."
    conda: ENVDIR + "/IMP_mapping.yaml"
    shell:
        """
        samtools view -b -F 5 -L {input[1]} {input[0]} 2>> {log} | samtools bam2fq -n - > {output} 2>> {log} 
        samtools view -b -f 65 -F 4 -L {input[1]} {input[0]} 2>> {log} | samtools bam2fq -n - >> {output} 2>> {log}
        samtools view -b -f 129 -F 4 -L {input[1]} {input[0]} 2>> {log} | samtools bam2fq -n - >> {output} 2>> {log}
        """
