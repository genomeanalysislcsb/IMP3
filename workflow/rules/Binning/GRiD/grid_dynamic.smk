
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
            report("Stats/%s/%s.bins.tsv" %(ASS,ASS),"Binning")
        resources:
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
            dir="Binning/MetaWrap/bins",
            ctrl="Stats/%s/%s.bins.tsv" %(ASS,ASS),
            bins=gather_bins
        output:
            "Binning/per_bin_results.tar.gz"
        threads: 1
        resources:
            runtime = "8:00:00",
            mem = MEMCORE
        log: "logs/binning_bins_tar.log"
        message: "bins_tar: Compressing per-bin results of GRiD."
        shell:
            """
            if [ -z "$(ls -A {input.dir})" ]; then
              tar cvzf {output} {input.bins} >> {log} 2>&1 && rm -r {input.bins} >> {log} 2>&1
            else
              touch {output}
            fi
            """



checkpoint getBins:
    input: 
        "Binning/selected_DASTool_summary.txt"
    output:
        "status/binning_GRiD.prepped"
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
#    log: "logs/binning_prepareGRiD.log"
    message: "getBins: Copying MetaWrap bins to single directories."
    run:
        if os.stat("Binning/MetaWrap/metawrap_bins.stats").st_size > 0:
            selected_bins=pd.read_table("Binning/MetaWrap/metawrap_bins.stats")['bin']
            for BIN in selected_bins:

                shell("mkdir -p Binning/MetaWrap/bins/{BIN} && cp Binning/MetaWrap/bins/{BIN}.fasta Binning/MetaWrap/bins/{BIN}/contigs.fasta")
        shell("touch {output}")

rule grid:
    input:
       "Binning/MetaWrap/bins/{clusterID}/contigs.fasta",
       "Binning/MetaWrap/bins/{clusterID}/reads.fastq.gz"
    output:
        directory("Binning/MetaWrap/bins/{clusterID}/grid")
    params:
        clusterInput="Binning/MetaWrap/bins/"
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_grid.yaml"
    log: "logs/binning_grid.{clusterID}.log"
    message: "grid: Running GRiD on {wildcards.clusterID}."
    shell:
        """
        grid single -r {params.clusterInput}{wildcards.clusterID} -e fastq.gz -o {output} -g {input[0]} > {log} 2>&1
        """

rule prepare_gridBed:
    input:
        "Binning/MetaWrap/bins/{clusterID}/contigs.fasta"
    output:
        "Binning/MetaWrap/bins/{clusterID}/contigs.bed"
    resources:
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
        "Binning/MetaWrap/bins/{clusterID}/contigs.bed"
    output:
        "Binning/MetaWrap/bins/{clusterID}/reads.fastq.gz"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: 1
    log: "logs/binning_prepare_gridReads.{clusterID}.log"
    message: "prepare_gridReads: Preparing fastq file for GRiD - {wildcards.clusterID}."
    conda: ENVDIR + "/IMP_mapping.yaml"
    shell:
        """
        samtools view -b -F 5 -L {input[1]} {input[0]} 2>> {log} | samtools fastq -n - > {output} 2>> {log} 
        samtools view -b -f 65 -F 4 -L {input[1]} {input[0]} 2>> {log} | samtools fastq -n - >> {output} 2>> {log}
        samtools view -b -f 129 -F 4 -L {input[1]} {input[0]} 2>> {log} | samtools fastq -n - >> {output} 2>> {log}
        """
