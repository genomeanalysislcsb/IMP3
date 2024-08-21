# TODO: Manually changed comebin cluster.py to try to deal with cases where no seed clusters are found

rule comebin_setup_scripts_from_fork:
    output:
        comebin_path = directory(os.path.join(ROOTDIR, "submodules/comebin")), 
        comebin_ready_check = os.path.join(ROOTDIR, 'submodules/comebin.ready')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        4
    conda:
        os.path.join(ENVDIR, "IMP_comebin.yaml")
    log:
        os.path.join(ENVDIR, "logs/comebin_setup_scripts.log")
    params:
        fork_url = "https://github.com/ohickl/COMEBin.git",
        fork_branch = "improve_cluster_leiden"
    message:
        "comebin: Setting up scripts from fork."
    shell:
        """
        # Clone git repo to submodule dir
        echo "Cloning comebin fork ..." >> {log}
        git clone -b {params.fork_branch} {params.fork_url} {output.comebin_path} >> {log} 2>&1

        # Replace all instance of main.py, cluster.py and utils.py with the ones from the fork
        # echo "Replacing scripts in conda env ..." >> {log}
        # find ${{CONDA_PREFIX}}/bin/COMEBin -type f -name main.py -exec cp {output.comebin_path}/COMEBin/main.py {{}} \; >> {log} 2>&1
        find ${{CONDA_PREFIX}}/bin/COMEBin -type f -name cluster.py -exec cp {output.comebin_path}/COMEBin/cluster.py {{}} \; >> {log} 2>&1
        find ${{CONDA_PREFIX}}/bin/COMEBin -type f -name utils.py -exec cp {output.comebin_path}/COMEBin/utils.py {{}} \; >> {log} 2>&1

        # Copy new run scripts into conda env bin as well as conda env comebin dir
        echo "Copying run scripts to conda env bin ..." >> {log}
        cp {output.comebin_path}/bin/run_comebin_train.sh ${{CONDA_PREFIX}}/bin/ >> {log} 2>&1
        cp {output.comebin_path}/bin/run_comebin_bin.sh ${{CONDA_PREFIX}}/bin/ >> {log} 2>&1

        echo "Copying run scripts to conda env comebin dir ..." >> {log}
        # Find comebin conda env dir
        comebin_env_dir=$(find ${{CONDA_PREFIX}} -type d -name COMEBin)
        cp {output.comebin_path}/bin/run_comebin_train.sh ${{comebin_env_dir}}/ >> {log} 2>&1
        cp {output.comebin_path}/bin/run_comebin_bin.sh ${{comebin_env_dir}}/ >> {log} 2>&1

        # Finally, we need to resintall torch using pip, to avoid unhinged library linkage problems for now
        pip3 install torch --force-reinstall >> {log} 2>&1

        # Touch ready file
        touch {output.comebin_ready_check} >> {log} 2>&1

        echo "Done." >> {log}
        """

def set_runtime_train(wildcards, attempt):
    return 480 if attempt == 1 else LONG_RUNTIME

rule comebin_train:
    input:
        comebin_ready_check = os.path.join(ROOTDIR, 'submodules/comebin.ready'),
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        mapping = "Assembly/mg.reads.sorted.bam"
    output:
        out_dir = directory("Binning/COMEBin"),
        filt_ass = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa'),
        filt_ass_lst = temp(f'Assembly/contig_length_filter1000.txt'),
        filt_ass_len = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa_lengths.txt'),
        tmp_align_dir = temp(directory('Binning/COMEBin/bam')),
        tmp_feature_dir = temp(directory('Binning/COMEBin/data_augmentation')),
        tmp_res_dir = temp(directory('Binning/COMEBin/comebin_res')),
        embeddings = temp('Binning/COMEBin/comebin_res/embeddings.tsv'),
        ass_stats = temp("Binning/COMEBin/assembly_stats.tsv")
    params:
        abs_out_dir = OUTPUTDIR,
        min_cont_size = config["binning"]["comebin"]["cutoff"],
        impossible_check = 'status/binning_comebin.impossible'
    threads:
        getThreads(7, high_mem_cluster=True)
    resources:
        runtime = set_runtime_train,
        mem_mb = 27000
    retries:
        2
    conda:
        os.path.join(ENVDIR, "IMP_comebin.yaml")
    log:
        os.path.join(OUTPUTDIR, "logs/binning/COMEBin_train.log")
    message:
        "Running COMEBin training."
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
          touch {output.filt_ass} \
                {output.filt_ass_lst} \
                {output.filt_ass_len} \
                {output.embeddings} \
                {params.impossible_check} >> {log} 2>&1
          mkdir -p {output.out_dir} \
                   {output.tmp_align_dir} \
                   {output.tmp_feature_dir} \
                   {output.tmp_res_dir} >> {log} 2>&1
          exit 0
        fi
        
        # Change to outdir to avoid weird errors
        echo "Moving to {params.abs_out_dir}/{output.out_dir} ..." >> {log} 2>&1 
        cd {params.abs_out_dir}/{output.out_dir} >> {log} 2>&1

        # Prepare tmp dir for mg read alignments
        echo "Symlinking mg read alignments to tmp dir {params.abs_out_dir}/{output.tmp_align_dir} ..." >> {log} 2>&1
        mkdir -p {params.abs_out_dir}/{output.tmp_align_dir} >> {log} 2>&1
        ln -s {params.abs_out_dir}/{input.mapping} {params.abs_out_dir}/{output.tmp_align_dir}/ >> {log} 2>&1

        # Filter contigs by length
        echo "Filtering contigs by length ..." >> {log} 2>&1
        $CONDA_PREFIX/bin/COMEBin/scripts/Filter_tooshort.py {params.abs_out_dir}/{input.assembly} {params.min_cont_size} >> {log} 2>&1
        echo "Filtered assembly written to {params.abs_out_dir}/{output.filt_ass}." >> {log} 2>&1

        # Run COMEBin
        echo "Running COMEBin training ..." >> {log} 2>&1

        run_comebin_train.sh -a {params.abs_out_dir}/{output.filt_ass} \
                             -o {params.abs_out_dir}/{output.out_dir} \
                             -p {params.abs_out_dir}/{output.tmp_align_dir} \
                             -t {threads} >> {log} 2>&1

        echo "Done." >> {log} 2>&1
        """

# def set_runtime_bin(wildcards, attempt):
#     return LONG_RUNTIME if attempt == 1 else 7200

# def set_slurm_extra(wildcards, attempt):
#     return '--qos=normal' if attempt == 1 else '--qos=long'

rule comebin_bin:
    input:
        comebin_ready_check = os.path.join(ROOTDIR, 'submodules/comebin.ready'),
        out_dir = 'Binning/COMEBin',
        filt_ass = f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa',
        filt_ass_lst = f'Assembly/contig_length_filter1000.txt',
        filt_ass_len = f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa_lengths.txt',
        tmp_feature_dir = 'Binning/COMEBin/data_augmentation',
        tmp_res_dir = 'Binning/COMEBin/comebin_res',
        embeddings = 'Binning/COMEBin/comebin_res/embeddings.tsv',
    output:
        s2b = "Binning/COMEBin/contig2bin.tsv",
        out_bin_dir = directory("Binning/COMEBin/bins"),
        seed_file = temp(f'Assembly/{ASS}.assembly.merged_{config["binning"]["comebin"]["cutoff"]}.fa.bacar_marker.2quarter_lencutoff_1001.seed'),
    params:
        abs_out_dir = OUTPUTDIR,
        impossible_check = 'status/binning_comebin.impossible'
    threads:
        getThreads(128)
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB,
        # slurm_extra = set_slurm_extra
    conda:
        os.path.join(ENVDIR, "IMP_comebin.yaml")
    log:
        os.path.join(OUTPUTDIR, "logs/binning/COMEBin_bin.log")
    message:
        "Running COMEBin binning."
    shell:
        """
        # Check if assembly was large enough to train
        if [[ -f {params.impossible_check} ]]; then
          echo "### Less than 1000 contigs larger min size or assembly smaller than 2*10^6 bps. Exiting ... ###" >> {log}
          touch {output.s2b} \
                {output.seed_file} >> {log} 2>&1
          mkdir -p {output.out_bin_dir} >> {log} 2>&1
          exit 0
        fi

        # Change to outdir to avoid weird errors
        echo "Moving to {params.abs_out_dir}/{input.out_dir} ..." >> {log} 2>&1 
        cd {params.abs_out_dir}/{input.out_dir} >> {log} 2>&1

        # Run COMEBin binning
        echo "Running COMEBin binning ..." >> {log} 2>&1
        CUDA_VISIBLE_DEVICES=0

        run_comebin_bin.sh -o {params.abs_out_dir}/{input.out_dir} \
                           -a {params.abs_out_dir}/{input.filt_ass} \
                           -t {threads} >> {log} 2>&1

        # Move bin fastas and rename.
        echo "Moving and renaming bins ..." >> {log} 2>&1
        mv {params.abs_out_dir}/{input.out_dir}/comebin_res/comebin_res_bins {params.abs_out_dir}/{input.out_dir}/bins >> {log} 2>&1
        
        shopt -s nullglob
        for i in {params.abs_out_dir}/{output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {params.abs_out_dir}/{output.out_bin_dir}/COMEBIN_${{bin_name}}.fasta >> {log} 2>&1
        done
        shopt -u nullglob

        # Write bin table
        echo "Writing bin table ..." >> {log} 2>&1
        python {SRCDIR}/write_bin_table_from_fastas.py {params.abs_out_dir}/{output.out_bin_dir} {params.abs_out_dir}/{output.s2b} >> {log} 2>&1
        
        echo "Done." >> {log} 2>&1
        """
