rule metadecoder:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        bam = "Assembly/mg.reads.sorted.bam"
    output:
        out_dir = directory("Binning/MetaDecoder"),
        out_bin_dir = directory("Binning/MetaDecoder/bins"),
        ass_stats = temp("Binning/MetaDecoder/assembly_stats.tsv")
    params:
        min_cont_size = config["binning"]["MetaDecoder"]["cutoff"],
        impossible_check = f'{OUTPUTDIR}/status/binning_metadecoder.impossible'
    threads:
        getThreads(32)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metadecoder.yaml"
    log:
        os.path.join(OUTPUTDIR, "logs/binning/MetaDecoder.log")
    message:
        "Running MetaDecoder."
    shadow: "copy-minimal"
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
          touch {params.impossible_check} >> {log} 2>&1
          mkdir -p {output.out_dir} \
                   {output.out_bin_dir} >> {log} 2>&1
          exit 0
        fi

        # Get current work dir
        wdir=$(pwd -P)
        mkdir -p {output.out_dir} >> {log} 2>&1
        cd {output.out_dir} >> {log} 2>&1

        printf "MetaDecoder version: " >> {log} 2>&1
        metadecoder --version >> {log} 2>&1
        
        # Run MetaDecoder
        echo "Running MetaDecoder coverage ..." >> {log} 2>&1
        metadecoder coverage --bam $wdir/{input.bam} \
                             --threads {threads} \
                             --output $wdir/METADECODER.COVERAGE >> {log} 2>&1
        
        echo "Running MetaDecoder seed ..." >> {log} 2>&1
        metadecoder seed --fasta $wdir/{input.assembly} \
                         --threads {threads} \
                         --output $wdir/METADECODER.SEED >> {log} 2>&1

        echo "Running MetaDecoder cluster ..." >> {log} 2>&1
        metadecoder cluster --fasta $wdir/{input.assembly} \
                            --min_sequence_length {params.min_cont_size} \
                            --coverage $wdir/METADECODER.COVERAGE \
                            --seed $wdir/METADECODER.SEED \
                            --output $wdir/METADECODER >> {log} 2>&1

        # Move bin fastas and rename.
        echo "Moving bin fastas ..." >> {log} 2>&1
        mkdir -p $wdir/{output.out_bin_dir}
        
        shopt -s nullglob
        for i in $wdir/*.{{fa,fna,fasta}}; do
          bin_name=$(basename ${{i}} | sed 's/.[^.]*$//g')
          mv ${{i}} $wdir/{output.out_bin_dir}/${{bin_name}}.fasta >> {log} 2>&1
        done
        shopt -u nullglob

        echo "Done." >> {log} 2>&1
        """

metadecoder_bins = "Binning/MetaDecoder/bins"
metadecoder_s2b = "Binning/MetaDecoder/contig2bin.tsv"

rule metadecoder_contig2bin:
    input:
        out_bin_dir = "Binning/MetaDecoder/bins"
    output:
        metadecoder_s2b = "Binning/MetaDecoder/contig2bin.tsv"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "Preparing table for DAS_Tool from MetaDecoder output."
    run:
        import glob
        bin_fastas = glob.glob(f'{metadecoder_bins}/*.fasta')
        with open(metadecoder_s2b, 'a') as of:
            for bin_fasta in bin_fastas:
                bin_name = bin_fasta.split('/')[-1].replace('.fasta', '')
                with open(bin_fasta,'r') as bf:
                    for line in bf:
                        if line[0] == '>':
                            line_contig_id = line[1:].strip('\t \n')
                            of.write(f'{line_contig_id}\t{bin_name}\n')
