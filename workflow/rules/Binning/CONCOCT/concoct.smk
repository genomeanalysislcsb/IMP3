rule concoct:
    input:
        assembly="Assembly/%s.assembly.merged.fa" % ASS,
        depth="Stats/mg/%s.assembly.contig_depth.txt" % ASS
    output:
        out_dir = directory("Binning/CONCOCT"),
        s2b = "Binning/CONCOCT/contig2bin.tsv",
        out_bin_dir = directory("Binning/CONCOCT/bins"),
        ass_stats = temp("Binning/CONCOCT/assembly_stats.tsv")
    params:
        min_cont_size=config['binning']['CONCOCT']['cutoff'],
        impossible_check = 'status/binning_concoct.impossible'
    threads:
        getThreads(14)
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_concoct.yaml"
    log:
        "logs/binning/concoct.log"
    message:
        "concoct: Running CONCOCT."
    shell:
        """
        # Clear previous status check
        rm -rf {params.impossible_check} > {log} 2>&1

        # Check if there are enough contigs < min_cont_size to proceed with binning
        echo "### Checking number of contigs > min length ... ###" >> {log}
        # Stats format is `file  format  type  num_seqs    sum_len  min_len  avg_len  max_len`
        seqkit seq -m {params.min_cont_size} {input.assembly} 2>> {log} | seqkit stats -j 4 -T 2>> {log} | tail -n +2 > {output.ass_stats} 2>> {log}
        ass_n_contigs=$(cut -f 4 {output.ass_stats} 2>> {log})
        ass_size=$(cut -f 5 {output.ass_stats} 2>> {log})

        # Touch output and extit if theres less than 1000 contigs or the size is < 2*10^6 bps
        if [[ ! $ass_n_contigs -ge 1000 || ! $ass_size -ge 2000000 ]];then
          echo "### Less than 1000 contigs larger min size or assembly smaller than 2*10^6 bps. Exiting ... ###" >> {log}
          touch {output.s2b} \
                {params.impossible_check} >> {log} 2>&1
          mkdir -p {output.out_dir} \
                   {output.out_bin_dir} >> {log} 2>&1
          exit 0
        fi

        mkdir -p {output.out_bin_dir}
        concoct --threads {threads} \
                --composition_file {input.assembly} \
                --coverage_file <(cat <(printf "contig\tcov_mean_sample\n") {input.depth}) \
                -l {params.min_cont_size} \
                -b {output.out_dir}/ >> {log} 2>&1
                
        tree {output.out_dir} >> {log} 2>&1
                
        tail -n +2 <(sed -e 's/,/\t/g' {output.out_dir}/clustering_gt{params.min_cont_size}.csv) > {output.s2b}
        
        # Extract bins as individual FASTA
        extract_fasta_bins.py --output_path {output.out_bin_dir} {input.assembly} {output.out_dir}/clustering_gt{params.min_cont_size}.csv
        # Move bin fastas and rename.
        shopt -s nullglob
        for i in {output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
        done
        shopt -u nullglob
        rm {output.out_dir}/clustering_gt{params.min_cont_size}.csv
        """