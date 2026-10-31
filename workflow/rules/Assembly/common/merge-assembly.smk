if config["assembly"]["merge"] is None or config["assembly"]["merge"] == "":
    rule skip_merge_assembly:
        input:
            ass = 'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER)
        output:
            final_ass = 'Assembly/%s.assembly.merged.fa' % ASS,
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        params:
            sample = SAMPLE
        threads:
            6
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        log:
            f'logs/assembly/{ASS}.skip_merge.log'
        message:
            "skip_merge_assembly: Use the primary assembly only. Renaming contigs."
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
    rule merge_assembly_concat:
        input:
            ass_1 = 'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            ass_2 = 'Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            final_ass = 'Assembly/%s.assembly.merged.fa' % ASS,
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        params:
            sample = SAMPLE
        threads:
            6
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        log:
            f'logs/assembly/{ASS}.merge_concat.log'
        message:
            "merge_assembly_concat: Merging assemblies, renaming contigs."
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
elif config["assembly"]["merge"] == "cap3":
    rule merge_assembly_cap3:
        input:
            'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            'Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            'Assembly/%s.assembly.merged.fa' % ASS,
        resources:
            runtime = LONG_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        params:
            sample = SAMPLE
        threads: 1
        conda: os.path.join(ENVDIR, "IMP_assembly.yaml")
        log: "logs/assembly/merge_assembly_cap3.log"
        message: "merge_assembly_cap3: Merging assemblies."
        shadow: "copy-minimal"
        shell:
            """
            NAME_fin=Assembly/%s.assembly
            NAME=Assembly/intermediary/%s.assembly
            cat {input} > $NAME.cat.fa
            
            # Run cap3 
            cap3 $NAME.cat.fa -p {config[assembly][cap3][identity]} -o {config[assembly][cap3][overlap]} > {log} 2>&1
            
            # Concatenate assembled contigs, singletons and rename the contigs
            cat $NAME.cat.fa.cap.contigs $NAME.cat.fa.cap.singlets | \
            awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > $NAME_fin.merged.fa
            """ % (ASS,ASS)
elif config["assembly"]["merge"] == "flye":
    rule merge_assembly_flye:
        input:
            i1='Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            i2='Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            # temp files
            i1=temp("Assembly/intermediary/common_merge_flye/input.1.fasta"),
            i2=temp("Assembly/intermediary/common_merge_flye/input.2.fasta"),
            i12=temp("Assembly/intermediary/common_merge_flye/input.fasta"),
            uids=temp("Assembly/intermediary/common_merge_flye/unmapped.ids"),
            ufna=temp("Assembly/intermediary/common_merge_flye/unmapped.fasta"),
            # intermediate files
            flye="Assembly/intermediary/common_merge_flye/assembly.fasta",
            # final files
            final='Assembly/%s.assembly.merged.fa' % ASS,
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_BIGMEMCORE_MB
        params:
            sample = SAMPLE,
            iters = 0,
            minoverlap = 1000
        threads: getThreads(MAIN_BIGCORENO, bigmem=True)
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        log:
            "logs/assembly/merge_assembly_flye.log"
        message:
            "merge_assembly_flye: Merging assemblies."
        shadow: "copy-minimal"
        shell:
            """
            (
            # rename contigs in input files
            awk '/^>/{{print \">input1_contig_\" ++i; next}}{{print}}' {input.i1} > {output.i1}
            awk '/^>/{{print \">input2_contig_\" ++i; next}}{{print}}' {input.i2} > {output.i2}
            cat {output.i1} {output.i2} > {output.i12}
            # consensus assembly
            flye --subassemblies {output.i1} {output.i2} -o $(dirname {output.flye}) -t {threads} -i {params.iters} --min-overlap {params.minoverlap} 1>&2
            # align input contigs to assembly | get unmapped contigs | get contig IDs
            flye-minimap2 -ax asm5 {output.flye} {output.i12} | samtools view -f 4 - | cut -f1 > {output.uids}
            # filter input contigs by contig ID
            seqkit grep -w 60 -f {output.uids} {output.i12} -o {output.ufna}
            # add unmapped contigs | rename contigs
            cat {output.flye} {output.ufna} | awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > {output.final}
            ) &> {log}
            """
elif config["assembly"]["merge"] == "quickmerge":
    # TODO: Put merge_wrapper.py in the scripts folder. Document this
    rule format_fastas:
        input:
            i1 = 'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            i2 = 'Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER)
        output:
            if1 = "Assembly/intermediary/%s.%s_preprocessed.1_formatted.fa" % (ASS, IMP_ASSEMBLER),
            if2 = "Assembly/intermediary/%s.%s_unmapped.2_formatted.fa" % (ASS, IMP_ASSEMBLER)
        log:
            os.path.join(OUTPUTDIR,"logs","assembly", "merge_hybrid_format_fastas.log")
        benchmark:
            os.path.join(OUTPUTDIR,"logs","assembly", "merge_hybrid_format_fastas_benchmark.tsv")
        resources:
            runtime = VERY_SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE,
            # assembly_nr = lambda wildcards: wildcards.get("assembly_nr")
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
            assembly_1 = "Assembly/intermediary/%s.%s_preprocessed.1_formatted.fa" % (ASS, IMP_ASSEMBLER),
            assembly_2 = "Assembly/intermediary/%s.%s_unmapped.2_formatted.fa" % (ASS, IMP_ASSEMBLER)
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
            echo -e "File with contigs to join: {input.assembly_1}\\nFile with contigs to use as bridges: {input.assembly_2}." >> {log} 2>&1

            # Get current work dir
            wdir=$(pwd -P)

            mkdir -p {output.out_dir} >> {log} 2>&1
            cd {output.out_dir} >> {log} 2>&1

            {SRCDIR}/merge_wrapper.py \
                --threads {threads} \
                --version4 \
                --prefix {params.merge_prefix} \
                --merging_length_cutoff {params.min_overlap} \
                --length_cutoff {params.len_cutoff} \
                --hco {params.hco} \
                --c {params.c} \
                $wdir/{input.assembly_1} \
                $wdir/{input.assembly_2} >> {log} 2>&1

            parallel --jobs {threads} gzip ::: aln_summary_*.tsv \
                                            anchor_summary_*.txt \
                                            param_summary_*.txt \
                                            {params.merge_prefix}* >> {log} 2>&1

            seqkit replace -p .+ \
                        -r "S{params.sample}_{params.merge_prefix}_C{{nr}}" \
                        --nr-width 10 \
                        -o merged_{params.merge_prefix}_rn.fasta \
                        merged_{params.merge_prefix}.fasta >> {log} 2>&1

            rm -f merged_{params.merge_prefix}.fasta >> {log} 2>&1
            """

    use rule quickmerge_join_h1_with_h2 as quickmerge_join_h2_with_h1 with:
        input:
            assembly_1 = "Assembly/intermediary/%s.%s_unmapped.2_formatted.fa" % (ASS, IMP_ASSEMBLER),
            assembly_2 = "Assembly/intermediary/%s.%s_preprocessed.1_formatted.fa" % (ASS, IMP_ASSEMBLER)
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
            os.path.join(OUTPUTDIR, f"logs/assembly/{QM_MERGE_PREFIX_FINAL}.log")
        benchmark:
            os.path.join(OUTPUTDIR, f"logs/assembly/{QM_MERGE_PREFIX_FINAL}_benchmark.tsv")
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
            echo "File with contigs to join: {input.assembly_h1h2}." >> {log} 2>&1
            echo "File with contigs to use as bridges: {input.assembly_h2h1}." >> {log} 2>&1

            # Get current work dir
            wdir=$(pwd -P)

            mkdir -p {output.out_dir} >> {log} 2>&1

            cd {output.out_dir} >> {log} 2>&1

            {SRCDIR}/merge_wrapper.py \
                --threads {threads} \
                --version4 \
                --prefix {params.merge_prefix} \
                --merging_length_cutoff {params.min_overlap} \
                --length_cutoff {params.len_cutoff} \
                --hco {params.hco} \
                --c {params.c} \
                $wdir/{input.assembly_h1h2} \
                $wdir/{input.assembly_h2h1} >> {log} 2>&1

            parallel --jobs {threads} gzip ::: aln_summary_*.tsv \
                                            anchor_summary_*.txt \
                                            param_summary_*.txt \
                                            {params.merge_prefix}* >> {log} 2>&1

            # Rename final merges and "bridge" sequences
            seqkit replace -p .+ \
                        -r "S{params.sample}_merged_C{{nr}}" \
                        --nr-width 10 \
                        -o merged_{params.merge_prefix}_rn.fasta \
                        merged_{params.merge_prefix}.fasta >> {log} 2>&1

            seqkit replace -p .+ \
                        -r "S{params.sample}_bridges_C{{nr}}" \
                        --nr-width 10 \
                        -o bridges_h2h1_rn.fasta \
                        $wdir/{input.assembly_h2h1} >> {log} 2>&1

            cat merged_{params.merge_prefix}_rn.fasta \
                bridges_h2h1_rn.fasta \
                > merged_with_bridge_singletons.fasta 2>> {log}

            rm -f merged_{params.merge_prefix}_rn.fasta bridges_h2h1_rn.fasta >> {log} 2>&1

            echo "Running cd-hit-est." >> {log} 2>&1
            time cd-hit-est -i merged_with_bridge_singletons.fasta \
                            -o $wdir/{output.merged_assembly} \
                            -T {threads} \
                            -c 0.995 \
                            -d 50000 \
                            -M 300000 >> {log} 2>&1 \
            || time {params.alt_cdhit} -i merged_with_bridge_singletons.fasta \
                                    -o $wdir/{output.merged_assembly} \
                                    -T {threads} \
                                    -c 0.995 \
                                    -d 50000 \
                                    -M 300000 >> {log} 2>&1
            """
