import math

rule semibin:
    input:
        assembly = "Assembly/%s.assembly.merged.fa" % ASS,
        mapping = "Assembly/mg.reads.sorted.bam"
    output:
        s2b = "Binning/SemiBin/contig2bin.tsv",
        out_dir = directory("Binning/SemiBin"),
        out_bin_dir = directory("Binning/SemiBin/bins"),
        ass_stats = temp("Binning/SemiBin/assembly_stats.tsv")
    params:
        model = config["binning"]["SemiBin"]["model"],
        tmp_dir = os.path.join(config["tmp_dir"]),
        impossible_check = f'{OUTPUTDIR}/status/binning_semibin.impossible'
    threads:
        # get_max_cores_sb
        getThreads(128)
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB,
        cutoff = lambda wildcards, attempt: config["binning"]["SemiBin"]["cutoff"] + (250 * (attempt - 1))
    retries:
        5
    conda:
        ENVDIR + "/IMP_semibin.yaml"
    log:
        "logs/binning/SemiBin.log"
    message:
        "Running SemiBin2."
    shadow: "copy-minimal"
    shell:
        """
        # Clear previous status check
        rm -rf {params.impossible_check} > {log} 2>&1

        # Check if there are enough contigs < min_cont_size to proceed with binning
        echo "### Checking number of contigs > min length ... ###" >> {log}
        # Stats format is `file  format  type  num_seqs    sum_len  min_len  avg_len  max_len`
        seqkit seq -m {resources.cutoff} {input.assembly} 2>> {log} | seqkit stats -j 4 -T 2>> {log} | tail -n +2 > {output.ass_stats} 2>> {log}
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


        echo "Using cutoff: {resources.cutoff}"

        # Run SemiBin2
        SemiBin2 \
            single_easy_bin \
                --input-fasta {input.assembly} \
                --input-bam {input.mapping} \
                --environment {params.model} \
                --threads {threads} \
                --random-seed 0 \
                --min-len {resources.cutoff} \
                --compression none \
                --output {output.out_dir} >> {log} 2>&1

        # Move bin fastas and rename.
        mv {output.out_dir}/output_bins {output.out_dir}/bins
        shopt -s nullglob
        for i in {output.out_bin_dir}/*.fa; do
          bin_name=$(basename -s .fa ${{i}})
          mv ${{i}} {output.out_bin_dir}/${{bin_name}}.fasta
        done
        shopt -u nullglob

        # Write bin table
        python {SRCDIR}/write_bin_table_from_fastas.py {output.out_bin_dir} {output.s2b}
        rm -rf {output.out_dir}/output_bins
        """
