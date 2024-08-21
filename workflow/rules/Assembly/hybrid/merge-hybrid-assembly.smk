if config["assembly"]["merge"] is None or config["assembly"]["merge"] == "":
    localrules: skip_merge_hybrid_assembly
    rule skip_merge_hybrid_assembly:
        input:
            f'Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1.fa',
        output:
            f'Assembly/mgmt.assembly.merged.fa'
        params:
            sample = SAMPLE
        threads: 1
        message:
            "skip_merge_hybrid_assembly: Use the primary assembly only."
        shell:
            "awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {input} > {output}"
elif config["assembly"]["merge"] == "concat":
    localrules: merge_hybrid_assembly_concat
    rule merge_hybrid_assembly_concat:
        input:
            f'Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1.fa',
            f'Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2.fa',
        output:
            f'Assembly/mgmt.assembly.merged.fa'
        params:
            sample = SAMPLE
        threads: 1
        message:
            "merge_hybrid_assembly_concat: Merging assemblies."
        shell:
            "cat {input} | awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > {output}"
elif config["assembly"]["merge"] == "quickmerge":
    rule format_fastas:
        input:
            i1 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1.fa",
            i2 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2.fa"
        output:
            if1= f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1_formatted.fa",
            if2= f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2_formatted.fa"
        log:
            os.path.join(OUTPUTDIR,"logs","assembly", "merge_hybrid_format_fastas.log")
        benchmark:
            os.path.join(OUTPUTDIR,"logs","assembly", "merge_hybrid_format_fastas_benchmark.tsv")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE
        threads:
            4
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        message:
            "Formatting sample {params.sample} hybrid fastas."
        shell:
            """
            echo "Formatting {params.sample} hybrid fastas." >> {log}
            for i in {input}; do
                h_nr=${{i%.*}}
                h_nr=${{h_nr##*.}}
                seqkit replace -p .+ \
                               -r "S{params.sample}_H${{h_nr}}_C{{nr}}" \
                               --nr-width 10 \
                               -o ${{i/\.fa/_formatted.fa}} \
                               ${{i}} > {log} 2>&1
            done
            """

    # Need to find a better way ar at least move. Using params doesnt work
    QM_MERGE_PREFIX_FINAL = "h1h2_with_h2h1"
    QM_MERGE_PREFIX_1 = "h1_with_h2"
    QM_MERGE_PREFIX_2 = "h2_with_h1"

    rule quickmerge_join_h1_with_h2:
        input:
            qm_check = os.path.join(ENVDIR, "qm.ready"),
            assembly_1 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1_formatted.fa",
            assembly_2 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2_formatted.fa"
        output:
            merged_assembly = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_1}/merged_{QM_MERGE_PREFIX_1}_rn.fasta")
        log:
            os.path.join(OUTPUTDIR, "logs",  f"assembly/{QM_MERGE_PREFIX_1}.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  f"assembly/{QM_MERGE_PREFIX_1}_benchmark.tsv")
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            out_dir = directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_1}")),
            sample = SAMPLE,
            work_dir = OUTPUTDIR,
            min_overlap = config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff = config["assembly"]["quickmerge"]["len_cutoff"],
            hco = config["assembly"]["quickmerge"]["hco"],
            c = config["assembly"]["quickmerge"]["c"],
            merge_prefix = QM_MERGE_PREFIX_1
        threads:
            getThreads(128)
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        message:
            "Joining sample {params.sample} mgmt hybrid 1 with hybrid 2."
        shell:
            """
            echo -e "File with contigs to join: {input.assembly_1}.\nFile with contigs to use as bridges: {input.assembly_2}." >> {log} 2>&1
            mkdir -p {params.out_dir} >> {log} 2>&1

            # We have to change to the output directory because qickmerge will write the output files to the current directory ...
            cd {params.out_dir} >> {log} 2>&1

            # Run quickmerge
            time {SRCDIR}/merge_wrapper.py --threads {threads} \
                                  --version4 \
                                  --prefix {params.merge_prefix} \
                                  --merging_length_cutoff {params.min_overlap} \
                                  --length_cutoff {params.len_cutoff} \
                                  --hco {params.hco} \
                                  --c {params.c} \
                                  {params.work_dir}/{input.assembly_1} \
                                  {params.work_dir}/{input.assembly_2} >> {log} 2>&1

            # Archive the unused quickmerge output files
            parallel --jobs {threads} gzip ::: aln_summary_*.tsv anchor_summary_*.txt param_summary_*.txt {params.merge_prefix}* >> {log} 2>&1
            
            # Rename and format the merged assembly file
            seqkit replace -p .+ \
                           -r "S{params.sample}_{params.merge_prefix}_C{{nr}}" \
                           --nr-width 10 \
                           -o  {output.merged_assembly} \
                           {params.out_dir}/merged_{params.merge_prefix}.fasta >> {log} 2>&1
            
            # Remove the raw merged assembly file
            rm -f {params.out_dir}/merged_{params.merge_prefix}.fasta >> {log} 2>&1
            """

    use rule quickmerge_join_h1_with_h2 as quickmerge_join_h2_with_h1 with:
        input:
            assembly_1 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2_formatted.fa",
            assembly_2 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1_formatted.fa"
        output:
            merged_assembly = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_2}/merged_{QM_MERGE_PREFIX_2}_rn.fasta")
        log:
            os.path.join(OUTPUTDIR, "logs",  f"assembly/{QM_MERGE_PREFIX_2}.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  f"assembly/{QM_MERGE_PREFIX_2}_benchmark.tsv")
        params:
            out_dir = directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_2}")),
            sample = SAMPLE,
            work_dir = OUTPUTDIR,
            min_overlap = config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff = config["assembly"]["quickmerge"]["len_cutoff"],
            hco = config["assembly"]["quickmerge"]["hco"],
            c = config["assembly"]["quickmerge"]["c"],
            merge_prefix = QM_MERGE_PREFIX_2
        message:
            "Joining sample {params.sample} mgmt hybrid 2 with hybrid 1."

    rule quickmerge_join_h1h2_with_h2h1:
        input:
            assembly_h1h2 = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_1}/merged_{QM_MERGE_PREFIX_1}_rn.fasta"),
            assembly_h2h1 = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_2}/merged_{QM_MERGE_PREFIX_2}_rn.fasta")
        output:
            merged_assembly = f"Assembly/{ASS}.assembly.merged.fa"
        log:
            os.path.join(OUTPUTDIR, "logs",  f"assembly/{QM_MERGE_PREFIX_FINAL}.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  f"assembly/{QM_MERGE_PREFIX_FINAL}_benchmark.tsv")
        resources:
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            out_dir = directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", f"{QM_MERGE_PREFIX_FINAL}")),
            sample = SAMPLE,
            min_overlap = config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff = config["assembly"]["quickmerge"]["len_cutoff"],
            hco = config["assembly"]["quickmerge"]["hco"],
            c = config["assembly"]["quickmerge"]["c"],
            alt_cdhit = os.path.join(ROOTDIR, "submodules", "cdhit", "cd-hit-est"),
            merge_prefix = QM_MERGE_PREFIX_FINAL, 
            merge_prefix_1 = QM_MERGE_PREFIX_1,
            merge_prefix_2 = QM_MERGE_PREFIX_2
        threads:
            getThreads(128)
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        message:
            "Joining sample {params.sample} h1h2 merged contigs using h2h1 merged contigs."
        shell:
            """
            echo "File with contigs to join: {input.assembly_h1h2}." >> {log} 2>&1
            echo "File with contigs to use as bridges: {input.assembly_h2h1}." >> {log} 2>&1

            if [[ ! -d {params.out_dir} ]]; then
              mkdir -p {params.out_dir}
            fi

            cd {params.out_dir}

            time {SRCDIR}/merge_wrapper.py --threads {threads} \
                          --version4 \
                          --prefix {params.merge_prefix} \
                          --merging_length_cutoff {params.min_overlap} \
                          --length_cutoff {params.len_cutoff} \
                          --hco {params.hco} \
                          --c {params.c} \
                          {input.assembly_h1h2} \
                          {input.assembly_h2h1} >> {log} 2>&1
            
            for i in aln_summary_*.tsv anchor_summary_*.txt param_summary_*.txt {params.merge_prefix}*; do
              gzip ${{i}}
            done
            
            for fasta in {params.out_dir}/merged_{params.merge_prefix}.fasta {input.assembly_h2h1}; do
              if [[ $fasta == '{params.out_dir}/merged_{params.merge_prefix}.fasta' ]]; then
                out_file="{params.out_dir}/merged_{params.merge_prefix}_rn.fasta"
                suffix="merged"
              else
                out_file="{params.out_dir}/bridges_h2h1_rn.fasta"
                suffix="bridges"
              fi
              seqkit replace -p .+ \
                             -r "S{params.sample}_${{suffix}}_C{{nr}}" \
                             --nr-width 10 \
                             -o  ${{out_file}} \
                             ${{fasta}}
            done

            cat "{params.out_dir}/merged_{params.merge_prefix}_rn.fasta" \
                "{params.out_dir}/bridges_h2h1_rn.fasta" > "{params.out_dir}/merged_with_bridge_singletons.fasta"

            # cdhit
            echo "Running cdhit-est." >> {log} 2>&1
            time cd-hit-est -i {params.out_dir}/merged_with_bridge_singletons.fasta \
                            -o {OUTPUTDIR}/{output.merged_assembly} \
                            -T {threads} \
                            -c 0.995 \
                            -d 50000 \
                            -M 300000 >> {log} 2>&1 \
            || \
            time {params.alt_cdhit} -i {params.out_dir}/merged_with_bridge_singletons.fasta \
                                    -o {OUTPUTDIR}/{output.merged_assembly} \
                                    -T {threads} \
                                    -c 0.995 \
                                    -d 50000 \
                                    -M 300000 >> {log} 2>&1
            
            # rm -rf {params.out_dir}/merged_with_bridge_singletons.fasta {params.out_dir}/merged_{params.merge_prefix}.fasta \
            #        {params.out_dir}/bridges_h2h1_rn.fasta {params.out_dir}/merged_{params.merge_prefix}_rn.fasta
            """
