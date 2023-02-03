rule metadecoder:
    input:
        assembly="Assembly/%s.assembly.merged.fa" % ASS,
        bam="Assembly/mg.reads.sorted.bam"
    output:
        out_dir=directory("Binning/MetaDecoder"),
        out_bin_dir=directory("Binning/MetaDecoder/bins")
    params:
        min_cont_size=config["binning"]["MetaDecoder"]["cutoff"],
        work_dir = OUTPUTDIR
    threads:
        getThreads(8)
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    conda:
        ENVDIR + "/IMP_metadecoder.yaml"
    log:
        os.path.join(OUTPUTDIR, "logs/binning_MetaDecoder.log")
    message:
        "Running MetaDecoder."
    shell:
        """
        mkdir -p {params.work_dir}/{output.out_dir}
        cd {params.work_dir}/{output.out_dir}
        
        # Create sam files
        echo "Creating sam files." >> {log} 2>&1
        set -x
        metadecoder --version >> {log} 2>&1
        
        # parallel -j $THREADS "samtools view -@ 1 -h {{}} -o {params.work_dir}/{output.out_dir}/sam_files/{{/.}}.sam" ::: ${{ASSEMBLY_MAP[*]}}
        # SAM_FILES=({output.out_dir}/sam_files/*.sam)
        
        samtools view -@ {threads} -h {params.work_dir}/{input.bam} -o {params.work_dir}/{output.out_dir}/mg.reads.sorted.sam
        
        # Run MetaDecoder
        metadecoder coverage --sam {params.work_dir}/{output.out_dir}/mg.reads.sorted.sam \
                             --threads {threads} \
                             --output {params.work_dir}/{output.out_dir}/METADECODER.COVERAGE >> {log} 2>&1 \
        && \
        rm -rf {params.work_dir}/{output.out_dir}/mg.reads.sorted.sam \
        && \
        metadecoder seed --fasta {params.work_dir}/{input.assembly} \
                         --threads {threads} \
                         --output {params.work_dir}/{output.out_dir}/METADECODER.SEED >> {log} 2>&1 \
        && \
        metadecoder cluster --fasta {params.work_dir}/{input.assembly} \
                            --min_sequence_length {params.min_cont_size} \
                            --coverage {params.work_dir}/{output.out_dir}/METADECODER.COVERAGE \
                            --seed {params.work_dir}/{output.out_dir}/METADECODER.SEED \
                            --output {params.work_dir}/{output.out_dir}/METADECODER >> {log} 2>&1
        set +x
        
        # Move bin fastas and rename.
        mkdir {params.work_dir}/{output.out_dir}/bins
        mv {params.work_dir}/{output.out_dir}/*.fasta {params.work_dir}/{output.out_dir}/bins
        """

metadecoder_bins = "Binning/MetaDecoder/bins"
metadecoder_s2b = "Binning/MetaDecoder/scaffold2bin.tsv"

rule metadecoder_contig2bin:
    input:
        out_bin_dir="Binning/MetaDecoder/bins"
    output:
        metadecoder_s2b="Binning/MetaDecoder/scaffold2bin.tsv"
    resources:
        runtime = "00:10:00",
        mem_mb = MEMCORE_MB
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
