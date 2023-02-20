import math

rule GTDBtk_bins_single:
    input:
        bin = ancient("Binning/MetaWrap/bins/{clusterID}")
    output:
        gtdbtk_out = directory("Binning/MetaWrap/bins/{clusterID}/GTDB"),
        # summary_bac = "Binning/MetaWrap/bins/{clusterID}/GTDB/identify/gtdbtk.bac120.markers_summary.tsv",
        # summary_arc = "Binning/MetaWrap/bins/{clusterID}/GTDB/identify/gtdbtk.ar53.markers_summary.tsv",
    resources:
        runtime = "02:00:00",
        mem_mb = BIGMEMCORE_MB
    params:
        pplacer_cpus=getThreads(min(BIGCORENO, max(1, int(int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]] * BIGCORENO / 45))), bigmem=True),
        t1 = int(BIGMEMCORE[:-1]),
        t2 = to_gb_dict[config["mem"]["mem_unit"]],
        t3 = BIGCORENO,
        t4 = int(int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]] * BIGCORENO / 45),
    conda:
        os.path.join(ENVDIR, "IMP_gtdbtk.yaml")
    log:
        "logs/taxonomy_GTDBtk_bins_single.{clusterID}.log"
    message:
        "GTDBtk_bins_single: Classifiying {wildcards.clusterID} using GTDBtk."
    threads:
        getThreads(32, bigmem=True)
    shell:
        """
        export GTDBTK_DATA_PATH="{DBPATH}/gtdbtk/release207_v2"
        export PYTHONPATH=$CONDA_PREFIX/lib/python3.7/site-packages
        gtdbtk classify_wf --genome_dir {input.bin}/ \
                           --extension fasta \
                           --out_dir {output.gtdbtk_out} \
                           --cpus {threads} \
                           --pplacer_cpus {params.pplacer_cpus} >> {log} 2>&1 \
          || touch {output.gtdbtk_out}/gtdbtk.empty.summary.tsv
        """

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

rule GTDBtk_summarize:
    input:
        ancient(gather_GTDB_bins),
        mw_stats="Binning/MetaWrap/bins.stats"
    output:
        out_file=report("Stats/%s/%s.bins.tsv" % (ASS,ASS), category="Taxonomy"),
        grid="status/binning_GRiD.done" if MG else []
    resources:
        runtime = "00:15:00",
        mem_mb = MEMCORE_MB
    log:
        "logs/taxonomy_GTDBtk_summarize.log"
    message:
        "GTDBtk_summarize: Summarizing binning results."
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(10 / (int(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]]))), 1)
    run:
        # Check if MetaWrap output is empty
        if os.stat(input.mw_stats).st_size > 0:
            stats_data = {}

            # Get stats
            print('Getting MetaWrap bin stats.')
            with open(input.mw_stats,'r') as mws:
                # Get mw stats header
                header = next(mws).strip('\n \t')
                for line in mws:
                    line = line.strip('\n \t').split('\t')
                    bin = line[0]
                    stats_data[bin] = line

            # Read GRiD out
            got_grid_header = False

            print("Reading GRiD results.")
            # Get all GRiD stat files
            grid_out = glob.glob("Binning/MetaWrap/bins/*/grid/reads.GRiD.txt")
            if grid_out:
                for bin_grid_data in grid_out:
                    if os.stat(bin_grid_data).st_size > 0:
                        bin = bin_grid_data.split('/')[-3]
                        with open(bin_grid_data,'r') as f:
                            if not got_grid_header:
                                header = header + '\t' + '\t'.join(next(f).strip(' \t\n').split('\t')[1:])
                                got_grid_header = True
                            else:
                                next(f)
                            stats_data[bin].extend(next(f).strip(' \t\n').split('\t')[1:])
                if output.grid:
                    os.mknod(output.grid)

            # Read GTDBtk out
            got_gtdbtk_header = False

            print("Reading GTDBtk results.")
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
        else:
            print("No output from MetaWrap!")
            os.mknod(output.out_file)
            if output.grid:
                os.mknod(output.grid)

if MG:
    rule bins_tar:
        input:
            dir="Binning/MetaWrap/bins",
            ctrl="Stats/%s/%s.bins.tsv" % (ASS, ASS)
        output:
            "Binning/per_bin_results.tar.gz"
        threads:
            1
        resources:
            runtime = "8:00:00",
            mem_mb = MEMCORE_MB
        log:
            "logs/taxonomy_bins_tar.log"
        message:
            "bins_tar: Compressing per-bin results of GRiD and GTDBtk."
        shell:
            """
            if [ -z "$(ls -A {input.dir})" ]; then
              tar czf {output} {input.dir}/*/grid {input.dir}/*/GTDB && rm -r {input.dir}/*/grid {input.dir}/*/GTDB >> {log} 2>&1
            else
              touch {output}
            fi
            """
else:
    rule bins_tar:
        input:
            dir="Binning/MetaWrap/bins",
            ctrl="Stats/%s/%s.bins.tsv" % (ASS,ASS)
        output:
            "Binning/per_bin_results.tar.gz"
        threads:
            1
        resources:
            runtime = "8:00:00",
            mem_mb = MEMCORE_MB
        log:
            "logs/taxonomy_bins_tar.log"
        message:
            "bins_tar: Compressing per-bin results of GTDBtk."
        shell:
            """
            if [ -z "$(ls -A {input.dir})" ]; then
              tar czf {output} {input.dir}/*/GTDB && rm -r {input.dir}/*/GTDB >> {log} 2>&1
            else
              touch {output}
            fi
            """

