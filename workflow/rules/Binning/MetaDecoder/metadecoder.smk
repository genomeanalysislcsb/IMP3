rule metadecoder:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        bam = "Assembly/mg.reads.sorted.bam"
    output:
        out_dir = directory("Binning/MetaDecoder"),
        out_bin_dir = directory("Binning/MetaDecoder/bins")
    params:
        min_cont_size = config["binning"]["MetaDecoder"]["cutoff"],
        work_dir = OUTPUTDIR
    threads:
        getThreads(8)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metadecoder.yaml"
    log:
        os.path.join(OUTPUTDIR, "logs/binning/MetaDecoder.log")
    message:
        "Running MetaDecoder."
    shell:
        """
        mkdir -p {params.work_dir}/{output.out_dir}
        cd {params.work_dir}/{output.out_dir}

        printf "MetaDecoder version: " >> {log} 2>&1
        metadecoder --version >> {log} 2>&1

        # Create sam files
        echo "Creating sam files ..." >> {log} 2>&1
        samtools view -@ {threads} -h {params.work_dir}/{input.bam} -o {params.work_dir}/{output.out_dir}/mg.reads.sorted.sam >> {log} 2>&1
        
        # Run MetaDecoder
        echo "Running MetaDecoder coverage ..." >> {log} 2>&1
        metadecoder coverage --sam {params.work_dir}/{output.out_dir}/mg.reads.sorted.sam \
                             --threads {threads} \
                             --output {params.work_dir}/{output.out_dir}/METADECODER.COVERAGE >> {log} 2>&1
        
        rm -rf {params.work_dir}/{output.out_dir}/mg.reads.sorted.sam \
        
        echo "Running MetaDecoder seed ..." >> {log} 2>&1
        metadecoder seed --fasta {params.work_dir}/{input.assembly} \
                         --threads {threads} \
                         --output {params.work_dir}/{output.out_dir}/METADECODER.SEED >> {log} 2>&1

        echo "Running MetaDecoder cluster ..." >> {log} 2>&1
        metadecoder cluster --fasta {params.work_dir}/{input.assembly} \
                            --min_sequence_length {params.min_cont_size} \
                            --coverage {params.work_dir}/{output.out_dir}/METADECODER.COVERAGE \
                            --seed {params.work_dir}/{output.out_dir}/METADECODER.SEED \
                            --output {params.work_dir}/{output.out_dir}/METADECODER >> {log} 2>&1
        
        # Move bin fastas and rename.
        echo "Moving bin fastas ..." >> {log} 2>&1
        mkdir -p {params.work_dir}/{output.out_bin_dir}
        
        shopt -s nullglob
        for i in {output.out_dir}/*.{{fa,fna,fasta}}; do
          bin_file=${{i##*/}}
          bin_name=${{i%.*}}
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
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
