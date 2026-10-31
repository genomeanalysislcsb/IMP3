
checkpoint DAS_tool_check:
    input:
        expand("Binning/{binner}/contig2bin.tsv",binner=BIN_STEPS)
    output:
        "status/DAS_tool.ready"
    threads:
        getThreads(1)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    shadow: "minimal"
    shell:
        "touch {output}"

rule DAS_tool:
    input:
        contigs="Assembly/%s.assembly.merged.fa" % ASS,
        proteins="Binning/annotation.renamed.faa", 
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
        out_name="Binning/DASTool/selected",
        impossible_check = f'{OUTPUTDIR}/status/binning_DAStool.impossible'
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    conda:
        ENVDIR + "/IMP_annotation.yaml"
    log:
        "logs/binning/DAS_tool.log"
    message:
        "DAS_tool: Running DAS_Tool."
    shadow: "copy-minimal"
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
          touch {output.dt_summary} {output.s2b} {params.impossible_check}
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
          touch {output.dt_summary} {output.s2b} {params.impossible_check}
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
          touch {output.dt_summary} {output.s2b} {params.impossible_check}
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
        proteins = "Analysis/annotation/bakta/annotation.faa",
        annotations = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        contig2id = "Analysis/annotation/bakta/annotation.filt.contig2ID.tsv",
        tmp_proteins = temp("Binning/annotation.tmp.faa"),
        renamed_proteins = "Binning/annotation.renamed.faa"
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
    shadow: "copy-minimal"
    shell:
        """
        echo "Getting contig to ID mapping ..." > {log}
        cut -f 1,9 {input.annotations} 2>> {log} \
          | sort -k 1 2>> {log} \
          | perl -F'\\t' -ane '$_=~/ID=([^_]+)_([^;]+);/; print "$F[0]\\t$1\\t$2\\n";' > {output.contig2id} 2>> {log}

        # Exit with error if no contig2id file is created or it is empty
        if [[ ! -s {output.contig2id} ]]; then
          echo "No contig2id file created. Exiting." >> {log}
          exit 1
        fi

        echo "Renaming protein sequences ..." >> {log}
        sed 's/ .*//' {input.proteins} >> {output.tmp_proteins} 2>> {log}

        echo "Running prokka2dastoolHeader.pl ..." >> {log}
        {SRCDIR}/prokka2dastoolHeader.pl {output.tmp_proteins} \
                                         {output.contig2id} > {output.renamed_proteins} 2>> {log}

        # Exit with error if no renamed protein file is created or it is empty
        if [[ ! -s {output.renamed_proteins} ]]; then
          echo "No renamed protein file created. Exiting." >> {log}
          exit 1
        fi

        echo 'Done.'
        """

