rule concoct:
    input:
       assembly="Assembly/%s.assembly.merged.fa" % ASS,
       depth="Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        out_dir=directory("Binning/CONCOCT"),
        s2b="Binning/CONCOCT/scaffold2bin.tsv",
        out_bin_dir=directory("Binning/CONCOCT/bins")
    params:
        min_cont_size=config['binning']['CONCOCT']['cutoff']
    threads:
        getThreads(8)
    resources:
        runtime = "24:00:00",
        mem = MEMCORE
    conda:
        ENVDIR + "/IMP_metawrap.yaml"
    log:
        "logs/binning_concoct.log"
    message:
        "concoct: Running CONCOCT."
    shell:
        """
        mkdir -p {output.out_bin_dir}
        concoct --threads {threads} \
                --composition_file {input.assembly} \
                --coverage_file <(cat <(printf "contig\tcov_mean_sample\n") {input.depth}) \
                -l {params.min_cont_size} \
                -b {output.out_dir}/ >> {log} 2>&1
        tail -n +2 <(sed -e 's/,/\t/g' {output.out_dir}/clustering_gt1000.csv) > {output.s2b}
        
        # Extract bins as individual FASTA
        extract_fasta_bins.py --output_path {output.out_bin_dir} {input.assembly} {output.out_dir}/clustering_gt1000.csv
        # Move bin fastas and rename.
        shopt -s nullglob
        for i in {output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
        done
        shopt -u nullglob
        rm {output.out_dir}/clustering_gt1000.csv
        """