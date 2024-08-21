import glob
import math
import os

if MG and "taxonomy" not in IMP_STEPS:

    localrules: grid_all

    rule grid_all:
        input:
            gather_bins,
            # directory("Binning/MetaWrap/bins/{clusterID}/grid")
        output:
            "status/binning_GRiD.done"
        shell:
            """
            touch {output}
            """

    rule GRID_summarize:
        input:
            ctrl="status/binning_GRiD.done",
            mw_stats="Binning/MetaWrap/bins.stats"
        output:
            out_file=report("Stats/%s/%s.bins.tsv" %(ASS,ASS), "Binning")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        log:
            "logs/binning/GRID_summarize.log"
        message:
            "GRID_summarize: Summarizing binning results."
        threads:
            # Ensure at least 5gb of mem
            max(int(math.ceil(20 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
        run:
            # Check if MetaWrap output is empty
            if os.stat(input.mw_stats).st_size > 0:
                mw_stats_data = {}
                # Get stats
                print('Getting MetaWrap bin stats.')
                with open(input.mw_stats, 'r') as mws:
                    # Get mw stats header
                    header = next(mws).strip('\n \t')
                    for line in mws:
                        line = line.strip('\n \t').split('\t')
                        bin = line[0]
                        mw_stats_data[bin] = line

                # Read GRiD out
                grid_header_part_list = ['GRiD', '95% CI', 'GRiD unrefined', 'Species heterogeneity', 'Coverage', 'dnaA/ori ratio', 'ter/dif ratio']
                grid_header_part = '\t'.join(grid_header_part_list)

                print("Reading GRiD results.")
                with open(output.out_file, 'w') as of:
                    # Get all GRiD stat files
                    grid_out = glob.glob("Binning/MetaWrap/bins/*/grid")
                    if grid_out:
                        header = f'{header}\t{grid_header_part}\n'
                        of.write(header)
                        for bin_grid_data in grid_out:
                            bin = bin_grid_data.split('/')[-3]
                            grid_file = f'{grid_out}/reads.GRiD.txt'
                            if os.path.exists(grid_file) and not os.stat(grid_file).st_size == 0:
                                with open(grid_file, 'r') as f:
                                    # Skip header
                                    next(f)
                                    # Add and write single bin stat line
                                    of.write('\t'.join(mw_stats_data[bin] + next(f).split('\t')[1:]))
                            else:
                                # If no grid out for bin, write n tabs, depending on number of grid fields
                                of.write('\t'.join(mw_stats_data[bin] + ['\t'] * len(grid_header_part_list)))
                    else:
                        of.write(header)
                        for bin in mw_stats_data.keys():
                            of.write('\t'.join(mw_stats_data[bin] + ['\n']))
            else:
                print("No output from MetaWrap!")
                os.mknod(output)

    
    rule bins_tar:
        input:
            dir = ancient("Binning/MetaWrap/bins"),
            ctrl = ancient("Stats/%s/%s.bins.tsv" %(ASS,ASS))
            
        output:
            "Binning/per_bin_results.tar.gz"
        threads:
            getThreads(5)
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/misc/bins_tar.log"
        message:
            "bins_tar: Compressing per-bin results of GRiD."
        shell:
            """
            echo "Starting the compression process." >> {log}
            if find {input.dir} -type d -name grid -print -quit | grep -q .; then
                echo "Found grid directories to compress." >> {log}
                tar --use-compress-program=pigz -cvf {output} {input.dir}/*/grid >> {log} 2>&1
                if [ $? -eq 0 ]; then
                    echo "Compression successful. Now removing the original grid directories." >> {log}
                    find {input.dir} -type d -name grid -exec rm -r {} + >> {log} 2>&1
                    echo "Original grid directories removed successfully." >> {log}
                else
                    echo "Error during compression. Check the logs for details." >> {log}
                fi
            else
                echo "No grid directories found to compress. Creating an empty archive." >> {log}
                touch {output}
                echo "Empty archive created successfully." >> {log}
            fi
            echo "Compression process completed." >> {log}
            """

checkpoint getBins:
    input: 
        "Binning/MetaWrap/bins.stats"
    output:
        "status/binning_GRiD.prepped"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    log:
        "logs/binning/prepareGRiD.log"
    message:
        "getBins: Copying MetaWrap bins to single directories."
    shell:
        """
        # Check if bins.stats has a non-zero size
        if [ -s {input} ]; then
            # Read selected bins, copy them to single directories
            while read BIN; do
                if [ -f Binning/MetaWrap/bins/${{BIN}}.fasta ]; then
                    echo "Copying bin ${{BIN}} to single directory." >> {log} 2>&1
                    mkdir -p Binning/MetaWrap/bins/${{BIN}} >> {log} 2>&1 && \
                    cp Binning/MetaWrap/bins/${{BIN}}.fasta Binning/MetaWrap/bins/${{BIN}}/contigs.fasta >> {log} 2>&1
                fi
            done < <(tail -n +2 Binning/MetaWrap/bins.stats | cut -f1) >> {log} 2>&1
        fi

        # Touch the output file to signify completion
        touch {output} >> {log} 2>&1
        echo "Done." >> {log} 2>&1
        """

if MG:
    rule grid:
        input:
           checkp = "status/binning_GRiD.prepped",
           contigs="Binning/MetaWrap/bins/{clusterID}/contigs.fasta",
           reads="Binning/MetaWrap/bins/{clusterID}/reads.fastq.gz"
        output:
            directory("Binning/MetaWrap/bins/{clusterID}/grid")
        params:
            clusterInput="Binning/MetaWrap/bins"
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            # Ensure at least 5gb of mem
            max(int(math.ceil(20 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
        conda:
            ENVDIR + "/IMP_annotation.yaml"
        log:
            "logs/binning/grid/grid.{clusterID}.log"
        message:
            "grid: Running GRiD on {wildcards.clusterID}."
        shell:
            """
            grid single -r {params.clusterInput}/{wildcards.clusterID} \
                        -e fastq.gz \
                        -o {output} \
                        -g {input.contigs} > {log} 2>&1
            """

    rule prepare_gridBed:
        input:
            "Binning/MetaWrap/bins/{clusterID}/contigs.fasta"
        output:
            "Binning/MetaWrap/bins/{clusterID}/contigs.bed"
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            # Ensure at least 5gb of mem
            max(int(math.ceil(20 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
        log:
            "logs/binning/grid/prepare_gridBed.{clusterID}.log"
        message:
            "prepare_gridBed: Preparing .bed file for GRiD - {wildcards.clusterID}."
        shell:
            """
            {SRCDIR}/fastaBed.pl {input} >> {output} 2> {log}
            """

    rule prepare_gridReads:
        input:
            assembly="Assembly/mg.reads.sorted.bam",
            cont_cov="Binning/MetaWrap/bins/{clusterID}/contigs.bed"
        output:
            reads="Binning/MetaWrap/bins/{clusterID}/reads.fastq.gz"
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        threads:
            # Ensure at least 5gb of mem
            max(int(math.ceil(20 / (int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
        log:
            "logs/binning/grid/prepare_gridReads.{clusterID}.log"
        message:
            "prepare_gridReads: Preparing fastq file for GRiD - {wildcards.clusterID}."
        conda:
            ENVDIR + "/IMP_mapping.yaml"
        shell:
            """
            samtools view -b -F 5 -L {input.cont_cov} {input.assembly} 2>> {log} | \
              samtools fastq -n - > {output} 2>> {log}
            samtools view -b -f 65 -F 4 -L {input.cont_cov} {input.assembly} 2>> {log} | \
              samtools fastq -n - >> {output} 2>> {log}
            samtools view -b -f 129 -F 4 -L {input.cont_cov} {input.assembly} 2>> {log} | \
              samtools fastq -n - >> {output} 2>> {log}
            """
