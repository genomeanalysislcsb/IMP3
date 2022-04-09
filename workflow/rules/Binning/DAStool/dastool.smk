localrules: DAS_tool_check

checkpoint DAS_tool_check:
    input:
        expand("Binning/{binner}/scaffold2bin.tsv",binner=BIN_STEPS)
    output:
        "Binning/DAS_tool.ready"
    shell:
        "touch {output}"

rule DAS_tool:
    input:
        contigs="Assembly/%s.assembly.merged.fa" % ASS,
        proteins="Binning/prokka.renamed.faa", 
        tabs=check_bintabs
    output:
        dt_summary="Binning/DASTool/selected_DASTool_summary.txt",
        s2b="Binning/DASTool/selected_DASTool_scaffolds2bin.txt",
        bin_dir=directory("Binning/DASTool/bins"),
        out_dir=directory("Binning/DASTool")
    threads:
        getThreads(4)
    params:
        score_thresh=config["binning"]["DAS_Tool"]["score_threshhold"],
        out_name="Binning/DASTool/dastool"
    resources:
        runtime = "24:00:00",
        mem = BIGMEMCORE
    conda:
        ENVDIR + "/IMP_annotation.yaml"
    log:
        "logs/binning_DAS_tool.log"
    message:
        "DAS_tool: Running DAS_Tool."
    shell:
        """
        TABS="{input.tabs}"
        TABS_ARR=(${{TABS}})
        TAB_tmp=${{TABS//Binning\//}}
        TAB_tmp2=${{TAB_tmp//\/scaffold2bin.tsv/}}
        TABSl=${{TAB_tmp2// /,}}
        
        if [[ ! -d {output.out_dir} ]]; then
          echo "Creating DASTool out dir." >> {log} 2>&1
          mkdir -p {output.out_dir}
        fi
        
        if [ ${{#TABS_ARR[@]}} == 0 ]; then
          touch {output.dt_summary} {output.s2b} status/binning_DAStool.impossible
          mkdir {output.bin_dir}
          exit 0
        fi
        
        DAS_Tool --bins=${{TABS// /,}} \
                 --labels=${{TABSl}} \
                 --contigs={input.contigs} \
                 --search_engine=diamond \
                 --proteins={input.proteins} \
                 --write_bins \
                 --threads={threads} \
                 --score_threshold={params.score_thresh} \
                 --outputbasename={params.out_name} >> {log} 2>&1 ||\
                  touch {output.dt_summary} {output.s2b} status/binning_DAStool.impossible
    
        mv Binning/DASTool/dastool_DASTool_bins {output.bin_dir}

        # Rename bin fastas
        shopt -s nullglob
        for bin in {output.bin_dir}/*.fa; do
          bin_name=$(sed -r -e 's/.fasta|.fa//g' <(echo ${{bin}}))
          mv ${{bin}} ${{bin_name}}.fasta
        done
        shopt -u nullglob
        
        if [ ! -f {output.dt_summary} ]; then
          touch {output.dt_summary}  {output.s2b} status/binning_DAStool.impossible && mkdir -p  {output.bin_dir}
        fi
        """

rule format_protein_names:
    input:
        "Analysis/annotation/prokka.faa",
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Analysis/annotation/annotation.filt.contig2ID.tsv",
        "Binning/prokka.renamed.faa"
    resources:
        runtime = "4:00:00",
        mem = MEMCORE
    threads: 1
    log: "logs/binning_format_protein_names.log"
    message: "format_protein_names: Preparing protein sequences for DAS tool."
    shell:
        """
        TMP_PROT=$(mktemp --tmpdir={TMPDIR} -t "proteins_XXXXXX.faa")
        cut -f 1,9 {input[1]} | sort -k 1 | perl -F'\\t' -ane '$_=~/ID=(.+)_(.+);inference/; print "$F[0]\\t$1\\t$2\\n";' >{output[0]}
	    sed 's/ .*//' {input[0]} >> $TMP_PROT
        {SRCDIR}/prokka2dastoolHeader.pl  $TMP_PROT {output[0]} > {output[1]} 2> {log}
	"""

