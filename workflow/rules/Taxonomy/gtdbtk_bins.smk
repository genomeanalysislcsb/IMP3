import math
import logging

def setup_logger(logfile=None, level=logging.DEBUG):
    """
    Setup a logger
    """

    # Formatting
    formatter = logging.Formatter('%(asctime)s - %(message)s', '%d/%m/%Y %I:%M:%S %p')

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Console handler
    ch = logging.StreamHandler()
    logger.addHandler(ch)

    # File handler
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return

##########

localrules: GTDB_all

rule GTDB_all:
    input:
        "Binning/per_bin_results.tar.gz"
    output:
        "status/taxonomy_GTDBtk.done"
    shell:
        """
        touch {output}
        """

# checkpoint getBins:
#     input:
#         mw_stats="Binning/MetaWrap/bins.stats",
#         mw_bin_dir="Binning/MetaWrap/bins"
#     output:
#         out_check="status/binning_{clusterID}_GTDBtk.prepped",
#         bin_dirs=directory("Binning/MetaWrap/bins/{clusterID}"),
#         bin_dir_fastas="Binning/MetaWrap/bins/{clusterID}/contigs.fasta"
#     resources:
#         runtime=SHORT_RUNTIME,
#         mem_mb=MAIN_MEMCORE_MB
#     threads: 1
#     log:
#         "logs/binning/{clusterID}_prepare_GTDBtk.log"
#     message:
#         "getBins: Copying MetaWrap bin {wildcards.clusterID} to single directory."
#     shadow: "copy-minimal"
#     run:
#         import os
#         import shutil

#         # Paths for the bin
#         cluster_id = wildcards.clusterID
#         bin_dir = output.bin_dir
#         bin_fasta = output.bin_dir_fasta
#         original_fasta = os.path.join("Binning/MetaWrap/bins", f"{cluster_id}.fasta")

#         # Ensure the directory exists and copy the FASTA file
#         os.makedirs(bin_dir, exist_ok=True)
#         if os.path.isfile(original_fasta):
#             shutil.copy2(original_fasta, bin_fasta)
#         else:
#             raise FileNotFoundError(f"Missing FASTA file for bin: {original_fasta}")

#         # Touch the completion file
#         with open(output.out_check, "w") as status_file:
#             status_file.write("Completed")
# TODO: Get this right so we can actually specify and check for outputs
checkpoint getBins:
    input: 
        mw_stats = "Binning/MetaWrap/bins.stats"
    output:
        check = "status/binning_GTDBtk.prepped"
    params:
        imp_out_dir = OUTPUTDIR,
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    log:
        "logs/binning/prepare_GTDBtk.log"
    message:
        "getBins: Copying MetaWrap bins to single directories."
    shadow: "copy-minimal"
    shell:
        """
        # Check if bins.stats has a non-zero size
        if [ -s {input.mw_stats} ]; then
            echo "Preparing bin dirs." >> {log} 2>&1
            # Read selected bins, copy them to single directories
            while read BIN; do
                if [ -f {params.imp_out_dir}/Binning/MetaWrap/bins/${{BIN}}.fasta ]; then
                    echo "Copying bin ${{BIN}} to single directory." >> {log} 2>&1
                    mkdir -p {params.imp_out_dir}/Binning/MetaWrap/bins/${{BIN}} >> {log} 2>&1 && \
                    cp {params.imp_out_dir}/Binning/MetaWrap/bins/${{BIN}}.fasta \
                       {params.imp_out_dir}/Binning/MetaWrap/bins/${{BIN}}/contigs.fasta >> {log} 2>&1
                fi
            done < <(tail -n +2 {input.mw_stats} | cut -f1) >> {log} 2>&1
        else
            echo "MetaWrap stats file is empty." >> {log} 2>&1
        fi

        # Touch the output file to signify completion
        touch {output.check} >> {log} 2>&1
        echo "Done." >> {log} 2>&1
        """

rule prepare_bin_fast2bed:
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
        "logs/binning/fast2bed/{clusterID}.log"
    message:
        "prepare_gridBed: Preparing .bed file for stats summary - {wildcards.clusterID}."
    shadow: "copy-minimal"
    shell:
        """
        {SRCDIR}/fastaBed.pl {input} >> {output} 2> {log}
        """

rule GTDBtk_bins_single:
    input:
        bin = ancient("Binning/MetaWrap/bins/{clusterID}"),
        dir_prep = "status/binning_GTDBtk.prepped"
    output:
        gtdbtk_out = directory("Binning/MetaWrap/bins/{clusterID}/GTDB")
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    params:
        pplacer_cpus=getThreads(min(MAIN_BIGCORENO, max(1, int(int(MAIN_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]] * MAIN_BIGCORENO / 45))), bigmem=True),
        db_path = os.path.join(DBPATH, "gtdbtk/release220")
    conda:
        os.path.join(ENVDIR, "IMP_gtdbtk.yaml")
    log:
        "logs/taxonomy/GTDBtk/GTDBtk_bins_single.{clusterID}.log"
    message:
        "GTDBtk_bins_single: Classifiying {wildcards.clusterID} using GTDBtk."
    threads:
        getThreads(64, bigmem=True)
    shadow: "copy-minimal"
    shell:
        """
        export GTDBTK_DATA_PATH="{params.db_path}"
        export PYTHONPATH=${{CONDA_PREFIX}}/lib/python3.8/site-packages
        gtdbtk classify_wf --genome_dir {input.bin} \
                           --extension fasta \
                           --mash_db {params.db_path}/mash \
                           --mash_s 10000 \
                           --out_dir {output.gtdbtk_out} \
                           --cpus {threads} \
                           --pplacer_cpus {params.pplacer_cpus} >> {log} 2>&1 \
          || {{ mkdir -p {output.gtdbtk_out} && touch {output.gtdbtk_out}/gtdbtk.empty.summary.tsv ;}}
        """

rule GTDBtk_summarize:
    input:
        ancient(gather_GTDB_bins),
        mw_stats = "Binning/MetaWrap/bins.stats",
        dir_prep = "status/binning_GTDBtk.prepped",
    output:
        out_file = report("Stats/%s/%s.bins.tsv" % (ASS,ASS), category="Taxonomy"),
        check = "status/binning_GTDBtk.done" if MG else []
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    log:
        "logs/taxonomy/GTDBtk_summarize.log"
    message:
        "GTDBtk_summarize: Summarizing binning results."
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(10 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
    run:
        # Setup logger
        setup_logger(log[0], level=logging.INFO)

        # Check if MetaWrap output is empty
        if os.stat(input.mw_stats).st_size > 0:
            stats_data = {}

            # Get stats
            logging.info('Getting MetaWrap bin stats.')
            with open(input.mw_stats,'r') as mws:
                # Get mw stats header
                header = next(mws).strip('\n \t')
                for line in mws:
                    line = line.strip('\n \t').split('\t')
                    bin = line[0]
                    stats_data[bin] = line

            # # Read GRiD out
            # got_grid_header = False

            # logging.info("Reading GRiD results.")
            # # Get all GRiD stat files
            # grid_out = glob.glob("Binning/MetaWrap/bins/*/grid/reads.GRiD.txt")
            # if grid_out:
            #     for bin_grid_data in grid_out:
            #         if os.stat(bin_grid_data).st_size > 0:
            #             bin = bin_grid_data.split('/')[-3]
            #             with open(bin_grid_data,'r') as f:
            #                 if not got_grid_header:
            #                     header = header + '\t' + '\t'.join(next(f).strip(' \t\n').split('\t')[1:])
            #                     got_grid_header = True
            #                 else:
            #                     next(f)
            #                 stats_data[bin].extend(next(f).strip(' \t\n').split('\t')[1:])
            #     if output.grid:
            #         os.mknod(output.grid)

            # Read GTDBtk out
            got_gtdbtk_header = False

            logging.info("Reading GTDBtk results.")
            # Get all GTBtk stat files
            gtdbtk_out = glob.glob("Binning/MetaWrap/bins/*/GTDB/gtdbtk.*.summary.tsv")
            if gtdbtk_out:
                for bin_gtdbtk_data in gtdbtk_out:
                    if os.stat(bin_gtdbtk_data).st_size > 0:
                        bin = bin_gtdbtk_data.split('/')[-3]
                        with open(bin_gtdbtk_data,'r') as f:
                            if not got_gtdbtk_header:
                                header = header + '\t' + '\t'.join(next(f).strip(' \t\n').split('\t')[1:]) + '\n'
                                got_gtdbtk_header = True
                            else:
                                next(f)
                            stats_data[bin].extend(next(f).strip(' \t\n').split('\t')[1:])

            # Write output
            with open(output.out_file,'w') as of:
                of.write(header)
                for bin in stats_data.keys():
                    of.write('\t'.join(stats_data[bin]) + '\n')
            
            if output.check:
                os.mknod(output.check)
        else:
            logging.warning("No output from MetaWrap!")
            os.mknod(output.out_file)
            if output.check:
                os.mknod(output.check)

if MG:
    rule bins_tar:
        input:
            dir = ancient("Binning/MetaWrap/bins"),
            ctrl = ancient("Stats/%s/%s.bins.tsv" % (ASS, ASS))
        output:
            "Binning/per_bin_results.tar.gz"
        threads:
            getThreads(12)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/misc/bins_tar.log"
        message:
            "bins_tar: Compressing MG per-bin results of GTDBtk."
        shadow: "copy-minimal"
        shell:
            """
            echo "Starting the compression process for GRiD and/or GTDBtk results." >> {log}
            if find {input.dir} -type d \( -name grid -o -name GTDB \) -print -quit | grep -q .; then
                echo "Found grid and/or GTDB directories to compress." >> {log}
                # Use find to list directories and tar them safely, handling spaces and special characters
                find {input.dir} -type d \( -name grid -o -name GTDB \) -print0 | tar --use-compress-program=pigz -cvf {output} --null -T - >> {log} 2>&1
                if [ $? -eq 0 ]; then
                    echo "Compression successful. Now removing the original grid and/or GTDB directories." >> {log}
                    find {input.dir} -type d \( -name grid -o -name GTDB \) -exec rm -r "{{}}" + >> {log} 2>&1
                    echo "Original directories removed successfully." >> {log}
                else
                    echo "Error during compression. Check the logs for details." >> {log}
                fi
            else
                echo "No grid or GTDB directories found to compress. Creating an empty archive." >> {log}
                touch {output}
                echo "Empty archive created successfully." >> {log}
            fi
            echo "Compression process for GRiD and/or GTDBtk results completed." >> {log}
            """
else:
    rule bins_tar:
        input:
            bins_dir = ancient("Binning/MetaWrap/bins"),            
            ctrl = ancient("Stats/%s/%s.bins.tsv" % (ASS, ASS))
        output:
            bin_res_arc = "Binning/per_bin_results.tar.gz"
        threads:
            getThreads(12)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/misc/bins_tar.log"
        message:
            "bins_tar: Compressing MT per-bin results of GTDBtk."
        shadow: "copy-minimal"
        shell:
            """
            echo "Starting the compression process for GTDBtk results." >> {log}
            if find {input.bins_dir} -type d -name GTDB -print -quit | grep -q .; then
                echo "Found GTDB directories to compress." >> {log}
                tar --use-compress-program=pigz -cvf {output.bin_res_arc} {input.bins_dir}/*/GTDB >> {log} 2>&1
                if [ $? -eq 0 ]; then
                    echo "Compression successful. Now removing the original GTDB directories." >> {log}
                    find {input.bins_dir} -type d -name GTDB -exec rm -r {{}} + >> {log} 2>&1
                    echo "Original GTDB directories removed successfully." >> {log}
                else
                    echo "Error during compression. Check the logs for details." >> {log}
                fi
            else
                echo "No GTDB directories found to compress. Creating an empty archive." >> {log}
                touch {output.bin_res_arc}
                echo "Empty archive created successfully." >> {log}
            fi
            echo "Compression process for GTDBtk results completed." >> {log}
            """

