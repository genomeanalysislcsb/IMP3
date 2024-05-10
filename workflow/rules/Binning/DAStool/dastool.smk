localrules: DAS_tool_check

checkpoint DAS_tool_check:
    input:
        expand("Binning/{binner}/contig2bin.tsv",binner=BIN_STEPS)
    output:
        "Binning/DAS_tool.ready"
    shell:
        "touch {output}"

rule DAS_tool:
    input:
        contigs="Assembly/%s.assembly.merged.fa" % ASS,
        proteins="Binning/prokka.renamed.faa", 
        tabs=check_bintabs
        # tabs=expand("Binning/{binner}/contig2bin.tsv",binner=BIN_STEPS)  # if used without DAS_tool_check
    output:
        dt_summary = "Binning/DASTool/selected_DASTool_summary.txt",
        s2b = "Binning/DASTool/selected_DASTool_scaffolds2bin.txt",
        bin_dir = directory("Binning/DASTool/bins"),
        out_dir = directory("Binning/DASTool")
    threads:
        getThreads(4)
    params:
        score_thresh=config["binning"]["DAS_Tool"]["score_threshhold"],
        out_name="Binning/DASTool/selected"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    conda:
        ENVDIR + "/IMP_annotation.yaml"
    log:
        "logs/binning/DAS_tool.log"
    message:
        "DAS_tool: Running DAS_Tool."
    shell:
        """
        TABS="{input.tabs}"
        TABS_ARR=(${{TABS}})
        TAB_tmp=${{TABS//Binning\//}}
        TAB_tmp2=${{TAB_tmp//\/contig2bin.tsv/}}
        TABSl=${{TAB_tmp2// /,}}
        
        if [[ ! -d {output.out_dir} ]]; then
          echo "Creating DASTool out dir." >> {log} 2>&1
          mkdir -p {output.out_dir}
        fi
        
        ###
        # If the arry is empty or none of the `contig2bin.tsv` contain any bins, exit
        if [ ${{#TABS_ARR[@]}} == 0 ]; then
          echo "No bins to process. Exiting." >> {log} 2>&1
          touch {output.dt_summary} {output.s2b} status/binning_DAStool.impossible
          mkdir {output.bin_dir}
          exit 0
        fi

        all_empty=true
        # Iterate over each file in the array
        for file in "${{TABS_ARR[@]}}"; do
            # If the file is not empty, set the flag to false and break the loop
            if [ -s "${{file}}" ]; then
                all_empty=false
                break
            fi
        done

        # If all files are empty, exit
        if ${{all_empty}}; then
          echo "No bins to process. Exiting." >> {log} 2>&1
          touch {output.dt_summary} {output.s2b} status/binning_DAStool.impossible
          mkdir {output.bin_dir}
          exit 0
        fi
        ###

        echo "Running DAS_Tool in bins: ${{TABS_ARR[@]}}" >> {log} 2>&1

        DAS_Tool --bins ${{TABS// /,}} \
                 --labels ${{TABSl}} \
                 --contigs {input.contigs} \
                 --search_engine diamond \
                 --proteins {input.proteins} \
                 --write_bins 1 \
                 --threads {threads} \
                 --score_threshold {params.score_thresh} \
                 --debug \
                 --outputbasename {params.out_name} >> {log} 2>&1
        
        # If DAS_Tool completes but chooses no bins, create empty files and exit
        if [[ ! -f {output.dt_summary} ]]; then
          echo "DAS_Tool completed but chose no bins. Exiting." >> {log} 2>&1
          touch {output.dt_summary} {output.s2b} status/binning_DAStool.impossible
          mkdir {output.bin_dir}
          exit 0
        fi

        if [[ -d {output.out_dir}/selected_DASTool_bins ]]; then
          echo "Moving DASTool bins to bin dir." >> {log} 2>&1
          mv {output.out_dir}/selected_DASTool_bins {output.bin_dir}

          # Rename bin fastas
          shopt -s nullglob
          for bin in {output.bin_dir}/*.fa; do
            bin_name=$(sed -r -e 's/.fasta|.fa//g' <(echo ${{bin}}))
            mv ${{bin}} ${{bin_name}}.fasta
          done
          shopt -u nullglob
        else
          echo "No bins found." >> {log} 2>&1
          mkdir -p {output.bin_dir}
        fi
        echo "DAS_Tool completed." >> {log} 2>&1
        """

rule format_protein_names:
    input:
        "Analysis/annotation/prokka.faa",
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Analysis/annotation/annotation.filt.contig2ID.tsv",
        "Binning/prokka.renamed.faa"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        # Ensure at least 10gb of mem
        max(int(math.ceil(10 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))),1)
    params:
        tmp_dir=TMPDIR
    log:
        "logs/binning/format_protein_names.log"
    message:
        "format_protein_names: Preparing protein sequences for DAS tool."
    shell:
        """
        export TMPDIR="{params.tmp_dir}"
        TMP_PROT=$(mktemp "${{TMPDIR}}/proteins_XXXXXX.faa")
        cut -f 1,9 {input[1]} | sort -k 1 | perl -F'\\t' -ane '$_=~/ID=(.+)_(.+);inference/; print "$F[0]\\t$1\\t$2\\n";' >{output[0]}
        sed 's/ .*//' {input[0]} >> $TMP_PROT
        {SRCDIR}/prokka2dastoolHeader.pl  $TMP_PROT {output[0]} > {output[1]} 2> {log}
        """

