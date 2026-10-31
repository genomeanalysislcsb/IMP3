if config["assembly"]["merge"] is None or config["assembly"]["merge"] == "":
    rule skip_merge_hybrid_assembly:
        input:
            ass = f'Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1.fa',
        output:
            final_ass = f'Assembly/mgmt.assembly.merged.fa'
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        params:
            sample = SAMPLE
        threads:
            6
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        log:
            f'assembly/logs/mgmt.skip_merge.log'
        message:
            "skip_merge_hybrid_assembly: Use the primary assembly only. Renaming contigs."
        shadow: "copy-minimal"
        shell:
            """
            echo "Renaming contigs and writing to final file ..." > {log}
            seqkit replace -p .+ \
                        -r "S{params.sample}_C{{nr}}" \
                        --nr-width 10 \
                        -o {output.final_ass} \
                        {input.ass} >> {log} 2>&1
            echo "Done" >> {log}
            """
elif config["assembly"]["merge"] == "concat":
    rule merge_hybrid_assembly_concat:
        input:
            ass_1 = f'Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1.fa',
            ass_2 = f'Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2.fa',
        output:
            final_ass = f'Assembly/mgmt.assembly.merged.fa'
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        params:
            sample = SAMPLE
        threads:
            6
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        log:
            f'assembly/logs/mgmt.merge_concat.log'
        message:
            "merge_hybrid_assembly_concat: Merging assemblies, renaming contigs."
        shadow: "copy-minimal"
        shell:
            """
            echo "Concatenating assemblies, renaming contigs and writing to final file ..." > {log}
            cat {input.ass_1} {input.ass_2} 2>> {log} \
                | seqkit replace \
                -p .+ \
                -r "S{params.sample}_C{{nr}}" \
                --nr-width 10 \
                -o {output.final_ass} \
                - >> {log} 2>&1
            echo "Done" >> {log}
            """
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
        shadow: "copy-minimal"
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

    rule install_qm_mummer4:
        output:
            qm_check = os.path.join(ENVDIR, "qm.ready"),
        log:
            os.path.join(ENVDIR, "logs", "quickmerge_install.log")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        threads:
            1
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        message:
            "Fixing quickmerge wrapper."
        shadow: "copy-minimal"
        shell:
            """
            wget https://raw.githubusercontent.com/mahulchak/quickmerge/master/merge_wrapper.py \
                -O ${{CONDA_PREFIX}}/bin/merge_wrapper.py >> {log} 2>&1 \
            && chmod 750 ${{CONDA_PREFIX}}/bin/merge_wrapper.py \
            && touch {output}
            """

    rule quickmerge_join_h1_with_h2:
        input:
            qm_check = os.path.join(ENVDIR, "qm.ready"),
            assembly_1 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1_formatted.fa",
            assembly_2 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2_formatted.fa"
        output:
            out_dir = directory(f"Assembly/intermediary/{QM_MERGE_PREFIX_1}"),
            merged_assembly = f"Assembly/intermediary/{QM_MERGE_PREFIX_1}/merged_{QM_MERGE_PREFIX_1}_rn.fasta"
        log:
            os.path.join(OUTPUTDIR, f"logs/assembly/{QM_MERGE_PREFIX_1}.log")
        benchmark:
            os.path.join(OUTPUTDIR, f"logs/assembly/{QM_MERGE_PREFIX_1}_benchmark.tsv")
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE,
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
        shadow: "copy-minimal"
        shell:
            """
            # Get absolute path to input assemblies
            ass_h1=$(readlink -f {input.assembly_1})
            ass_h2=$(readlink -f {input.assembly_2})

            echo "File with contigs to join: $ass_h1." >> {log} 2>&1
            echo "File with contigs to use as bridges: $ass_h2." >> {log} 2>&1

            mkdir -p {output.out_dir} >> {log} 2>&1
            cd {output.out_dir} >> {log} 2>&1

            echo "Running quickmerge ..." >> {log} 2>&1
            {SRCDIR}/merge_wrapper.py \
                --threads {threads} \
                --version4 \
                --prefix {params.merge_prefix} \
                --merging_length_cutoff {params.min_overlap} \
                --length_cutoff {params.len_cutoff} \
                --hco {params.hco} \
                --c {params.c} \
                $ass_h1 \
                $ass_h2 >> {log} 2>&1

            echo "Compressing unused files ..." >> {log} 2>&1
            parallel --jobs {threads} gzip ::: aln_summary_*.tsv \
                                            anchor_summary_*.txt \
                                            param_summary_*.txt \
                                            {params.merge_prefix}* >> {log} 2>&1

            echo "Renaming contigs..." >> {log} 2>&1
            seqkit replace -p .+ \
                        -r "S{params.sample}_{params.merge_prefix}_C{{nr}}" \
                        --nr-width 10 \
                        -o merged_{params.merge_prefix}_rn.fasta \
                        merged_{params.merge_prefix}.fasta >> {log} 2>&1

            echo "Removing intermediary fasta ..." >> {log} 2>&1
            rm -f merged_{params.merge_prefix}.fasta >> {log} 2>&1

            echo "Done." >> {log} 2>&1
            """

    use rule quickmerge_join_h1_with_h2 as quickmerge_join_h2_with_h1 with:
        input:
            os.path.join(ENVDIR, "qm.ready"),
            assembly_1 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.2_formatted.fa",
            assembly_2 = f"Assembly/intermediary/mgmt.{IMP_ASSEMBLER}_hybrid.1_formatted.fa"
        output:
            out_dir = directory(f"Assembly/intermediary/{QM_MERGE_PREFIX_2}"),
            merged_assembly = f"Assembly/intermediary/{QM_MERGE_PREFIX_2}/merged_{QM_MERGE_PREFIX_2}_rn.fasta"
        log:
            os.path.join(OUTPUTDIR, f"logs/assembly/{QM_MERGE_PREFIX_2}.log")
        benchmark:
            os.path.join(OUTPUTDIR, f"logs/assembly/{QM_MERGE_PREFIX_2}_benchmark.tsv")
        params:
            sample = SAMPLE,
            min_overlap = config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff = config["assembly"]["quickmerge"]["len_cutoff"],
            hco = config["assembly"]["quickmerge"]["hco"],
            c = config["assembly"]["quickmerge"]["c"],
            merge_prefix = QM_MERGE_PREFIX_2
        message:
            "Joining sample {params.sample} mgmt hybrid 2 with hybrid 1."

    rule quickmerge_join_h1h2_with_h2h1:
        input:
            assembly_h1h2 = f"Assembly/intermediary/{QM_MERGE_PREFIX_1}/merged_{QM_MERGE_PREFIX_1}_rn.fasta",
            assembly_h2h1 = f"Assembly/intermediary/{QM_MERGE_PREFIX_2}/merged_{QM_MERGE_PREFIX_2}_rn.fasta"
        output:
            out_dir = directory(f"Assembly/intermediary/{QM_MERGE_PREFIX_FINAL}"),
            merged_assembly = f"Assembly/{ASS}.assembly.merged.fa"
        log:
            os.path.join(OUTPUTDIR,f"logs/assembly/{QM_MERGE_PREFIX_FINAL}.log")
        benchmark:
            os.path.join(OUTPUTDIR,f"logs/assembly/{QM_MERGE_PREFIX_FINAL}_benchmark.tsv")
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE,
            min_overlap = config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff = config["assembly"]["quickmerge"]["len_cutoff"],
            hco = config["assembly"]["quickmerge"]["hco"],
            c = config["assembly"]["quickmerge"]["c"],
            # Adjust cd-hit path relative to your workflow
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
        shadow: "copy-minimal"
        shell:
            """
            # Get absolute path to input assemblies
            work_dir=$(readlink -f .)
            ass_h1h2=$(readlink -f {input.assembly_h1h2})
            ass_h2h1=$(readlink -f {input.assembly_h2h1})

            echo "File with contigs to join: $ass_h1h2." >> {log} 2>&1
            echo "File with contigs to use as bridges: $ass_h2h1." >> {log} 2>&1

            mkdir -p {output.out_dir} >> {log} 2>&1
            cd {output.out_dir} >> {log} 2>&1

            echo "Running quickmerge ..." >> {log} 2>&1
            {SRCDIR}/merge_wrapper.py \
                --threads {threads} \
                --version4 \
                --prefix {params.merge_prefix} \
                --merging_length_cutoff {params.min_overlap} \
                --length_cutoff {params.len_cutoff} \
                --hco {params.hco} \
                --c {params.c} \
                $ass_h1h2 \
                $ass_h2h1 >> {log} 2>&1

            echo "Compressing unused files ..." >> {log} 2>&1
            parallel --jobs {threads} gzip ::: aln_summary_*.tsv \
                                            anchor_summary_*.txt \
                                            param_summary_*.txt \
                                            {params.merge_prefix}* >> {log} 2>&1

            # Rename final merged and "bridge" sequences
            echo "Renaming contigs..." >> {log} 2>&1
            seqkit replace -p .+ \
                        -r "S{params.sample}_merged_C{{nr}}" \
                        --nr-width 10 \
                        -o merged_{params.merge_prefix}_rn.fasta \
                        merged_{params.merge_prefix}.fasta >> {log} 2>&1

            seqkit replace -p .+ \
                        -r "S{params.sample}_bridges_C{{nr}}" \
                        --nr-width 10 \
                        -o bridges_h2h1_rn.fasta \
                        $ass_h2h1 >> {log} 2>&1

            echo "Concatenating merged with bridges ..." >> {log} 2>&1
            cat merged_{params.merge_prefix}_rn.fasta \
                bridges_h2h1_rn.fasta \
                > merged_with_bridge_singletons.fasta 2>> {log}

            echo "Removing intermediary fastas ..." >> {log} 2>&1
            rm -f merged_{params.merge_prefix}_rn.fasta bridges_h2h1_rn.fasta >> {log} 2>&1

            echo "Deduplicating concatenated fasta ..." >> {log} 2>&1
            cd-hit-est \
                -i merged_with_bridge_singletons.fasta \
                -o $work_dir/{output.merged_assembly} \
                -T {threads} \
                -c 0.995 \
                -d 50000 \
                -M 300000 >> {log} 2>&1 \
            || {params.alt_cdhit} \
                   -i merged_with_bridge_singletons.fasta \
                   -o $work_dir/{output.merged_assembly} \
                   -T {threads} \
                   -c 0.995 \
                   -d 50000 \
                   -M 300000 >> {log} 2>&1
        
            echo "Done." >> {log} 2>&1
            """
