if config["assembly"]["merge"] is None or config["assembly"]["merge"] == "":
    localrules: skip_merge_assembly
    rule skip_merge_assembly:
        input:
            'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER)
        output:
            'Assembly/%s.assembly.merged.fa' % ASS,
        params:
            sample = SAMPLE
        threads: 1
        message:
            "skip_merge_assembly: Use the primary assembly only."
        shell:
            "awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {input} > {output}"
elif config["assembly"]["merge"] == "concat":
    localrules: merge_assembly_concat
    rule merge_assembly_concat:
        input:
            'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            'Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            'Assembly/%s.assembly.merged.fa' % ASS,
        params:
            sample = SAMPLE
        threads: 1
        message:
            "merge_assembly_concat: Merging assemblies."
        shell:
            "cat {input} | awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > {output}"
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
    # Now already done during install
    # rule install_qm_mummer4:
    #     output:
    #         os.path.join(ENVDIR, "qm.ready")
    #     log:
    #         os.path.join(ENVDIR, "logs",  "assembly_quickmerge_install.log")
    #     benchmark:
    #         os.path.join(ENVDIR, "logs",  "assembly_quickmerge_install_benchmark.tsv")
    #     resources:
    #         runtime = VERY_SHORT_RUNTIME,
    #         mem_mb = MAIN_MEMCORE_MB
    #     threads:
    #         1
    #     conda:
    #         os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
    #     message:
    #         "Fixing quickmerge wrapper."
    #     shell:
    #         """
    #         wget https://raw.githubusercontent.com/mahulchak/quickmerge/master/merge_wrapper.py \
    #              -O ${{CONDA_PREFIX}}/bin/merge_wrapper.py >> {log} 2>&1 \
    #           && chmod 750 ${{CONDA_PREFIX}}/bin/merge_wrapper.py \
    #           && touch {output}
    #         """

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

    rule quickmerge_join_h1_with_h2:
        input:
            os.path.join(ENVDIR, "qm.ready"),
            assembly_1 = "Assembly/intermediary/%s.%s_preprocessed.1_formatted.fa" % (ASS, IMP_ASSEMBLER),
            assembly_2 = "Assembly/intermediary/%s.%s_unmapped.2_formatted.fa" % (ASS, IMP_ASSEMBLER)
        output:
            out_dir=directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h1_with_h2")),
            merged_assembly=os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h1_with_h2/merged_h1_with_h2_rn.fasta")
        log:
            os.path.join(OUTPUTDIR, "logs",  "assembly", "join_h1_with_h2.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  "assembly", "join_h1_with_h2_benchmark.tsv")
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE,
            work_dir = OUTPUTDIR,
            min_overlap=config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff=config["assembly"]["quickmerge"]["len_cutoff"],
            hco=config["assembly"]["quickmerge"]["hco"],
            c=config["assembly"]["quickmerge"]["c"]
        threads:
            getThreads(128)
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        message:
            "Joining sample {params.sample} %s hybrid 1 with hybrid 2." % ASS
        shell:
            """
            echo "File with contigs to join: {input.assembly_1}." >> {log} 2>&1
            echo "File with contigs to use as bridges: {input.assembly_2}." >> {log} 2>&1

            if [[ ! -d {output.out_dir} ]]; then
              mkdir -p {output.out_dir}
            fi

            cd {output.out_dir}
            time merge_wrapper.py --threads {threads} \
                          --version4 \
                          --prefix "h1_with_h2" \
                          --merging_length_cutoff {params.min_overlap} \
                          --length_cutoff {params.len_cutoff} \
                          --hco {params.hco} \
                          --c {params.c} \
                          {params.work_dir}/{input.assembly_1} \
                          {params.work_dir}/{input.assembly_2} >> {log} 2>&1 \
            && \
              for i in aln_summary_*.tsv anchor_summary_*.txt param_summary_*.txt h1_with_h2*; do
                gzip ${{i}}
              done \
            && \
              seqkit replace -p .+ \
                             -r "S{params.sample}_H1H2_C{{nr}}" \
                             --nr-width 10 \
                             -o  {output.merged_assembly} \
                             merged_h1_with_h2.fasta \
            && \
              rm -f merged_h1_with_h2.fasta
            """

    rule quickmerge_join_h2_with_h1:
        input:
            os.path.join(ENVDIR, "qm.ready"),
            assembly_1= "Assembly/intermediary/%s.%s_preprocessed.1_formatted.fa" % (ASS,IMP_ASSEMBLER),
            assembly_2= "Assembly/intermediary/%s.%s_unmapped.2_formatted.fa" % (ASS,IMP_ASSEMBLER)
        output:
            out_dir=directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h2_with_h1")),
            merged_assembly=os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h2_with_h1/merged_h2_with_h1_rn.fasta")
        log:
            os.path.join(OUTPUTDIR, "logs",  "assembly", "join_h2_with_h1.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  "assembly", "join_h2_with_h1_benchmark.tsv")
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE,
            work_dir = OUTPUTDIR,
            min_overlap=config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff=config["assembly"]["quickmerge"]["len_cutoff"],
            hco=config["assembly"]["quickmerge"]["hco"],
            c=config["assembly"]["quickmerge"]["c"]
        threads:
            getThreads(128)
        conda:
            os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
        message:
            "Joining sample {params.sample} %s hybrid 2 with hybrid 1." % ASS
        shell:
            """
            echo "File with contigs to join: {input.assembly_2}." >> {log} 2>&1
            echo "File with contigs to use as bridges: {input.assembly_1}." >> {log} 2>&1

            if [[ ! -d {output.out_dir} ]]; then
              mkdir -p {output.out_dir}
            fi

            cd {output.out_dir}
            time merge_wrapper.py --threads {threads} \
                          --version4 \
                          --prefix "h2_with_h1" \
                          --merging_length_cutoff {params.min_overlap} \
                          --length_cutoff {params.len_cutoff} \
                          --hco {params.hco} \
                          --c {params.c} \
                          {params.work_dir}/{input.assembly_2} \
                          {params.work_dir}/{input.assembly_1} >> {log} 2>&1 \
            && \
            for i in aln_summary_*.tsv anchor_summary_*.txt param_summary_*.txt h2_with_h1*; do
              gzip ${{i}}
            done \
            && \
            seqkit replace -p .+ \
                           -r "S{params.sample}_H2H1_C{{nr}}" \
                           --nr-width 10 \
                           -o  {output.merged_assembly} \
                           merged_h2_with_h1.fasta  \
            && \
            rm -f merged_h2_with_h1.fasta
            """

    rule quickmerge_join_h1h2_with_h2h1:
        input:
            os.path.join(ENVDIR, "qm.ready"),
            assembly_h1h2 = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h1_with_h2/merged_h1_with_h2_rn.fasta"),
            assembly_h2h1 = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h2_with_h1/merged_h2_with_h1_rn.fasta")
        output:
            out_dir=directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "h1h2_with_h2h1")),
            # merged_assembly=os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "h1h2_with_h2h1/merged_with_bridge_singletons_99id.fasta")
            merged_assembly = 'Assembly/%s.assembly.merged.fa' % ASS
        log:
            os.path.join(OUTPUTDIR, "logs",  "assembly", "h1h2_with_h2h1.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  "assembly", "h1h2_with_h2h1_benchmark.tsv")
        resources:
            runtime = MEDIUM_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        params:
            sample = SAMPLE,
            min_overlap=config["assembly"]["quickmerge"]["min_overlap"],
            len_cutoff=config["assembly"]["quickmerge"]["len_cutoff"],
            hco=config["assembly"]["quickmerge"]["hco"],
            c=config["assembly"]["quickmerge"]["c"],
            alt_cdhit = os.path.join(ROOTDIR, "submodules", "cdhit", "cd-hit-est")
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

            if [[ ! -d {output.out_dir} ]]; then
              mkdir -p {output.out_dir}
            fi

            cd {output.out_dir}

            time merge_wrapper.py --threads {threads} \
                          --version4 \
                          --prefix "h1h2_with_h2h1" \
                          --merging_length_cutoff {params.min_overlap} \
                          --length_cutoff {params.len_cutoff} \
                          --hco {params.hco} \
                          --c {params.c} \
                          {input.assembly_h1h2} \
                          {input.assembly_h2h1} >> {log} 2>&1 \
              && \
            for i in aln_summary_*.tsv anchor_summary_*.txt param_summary_*.txt h1h2_with_h2h1*; do
              gzip ${{i}}
            done \
              && \
            for fasta in {output.out_dir}/merged_h1h2_with_h2h1.fasta {input.assembly_h2h1}; do
              if [[ $fasta == '{output.out_dir}/merged_h1h2_with_h2h1.fasta' ]]; then
                out_file="{output.out_dir}/merged_h1h2_with_h2h1_rn.fasta"
                suffix="merged"
              else
                out_file="{output.out_dir}/bridges_h2h1_rn.fasta"
                suffix="bridges"
              fi
              seqkit replace -p .+ \
                             -r "S{params.sample}_${{suffix}}_C{{nr}}" \
                             --nr-width 10 \
                             -o  ${{out_file}} \
                             ${{fasta}}
            done

            cat "{output.out_dir}/merged_h1h2_with_h2h1_rn.fasta" \
                "{output.out_dir}/bridges_h2h1_rn.fasta" > "{output.out_dir}/merged_with_bridge_singletons.fasta"

            # cdhit
            echo "Running cdhit-est." >> {log} 2>&1
            time cd-hit-est -i {output.out_dir}/merged_with_bridge_singletons.fasta \
                            -o {OUTPUTDIR}/{output.merged_assembly} \
                            -T {threads} \
                            -c 0.995 \
                            -d 50000 \
                            -M 300000 >> {log} 2>&1 \
            || \
            time {params.alt_cdhit} -i {output.out_dir}/merged_with_bridge_singletons.fasta \
                                    -o {OUTPUTDIR}/{output.merged_assembly} \
                                    -T {threads} \
                                    -c 0.995 \
                                    -d 50000 \
                                    -M 300000 >> {log} 2>&1 \
            && \
            rm -rf {output.out_dir}/merged_with_bridge_singletons.fasta {output.out_dir}/merged_h1h2_with_h2h1.fasta \
                   {output.out_dir}/bridges_h2h1_rn.fasta {output.out_dir}/merged_h1h2_with_h2h1_rn.fasta \
            """

    # # Commented compression for now, as its messing with multi-step/re- running of imp
    # rule clean_up:
    #     input:
    #         ass1_form = "Assembly/intermediary/%s.%s_preprocessed.1_formatted.fa" % (ASS,IMP_ASSEMBLER),
    #         ass2_form = "Assembly/intermediary/%s.%s_unmapped.2_formatted.fa" % (ASS,IMP_ASSEMBLER),
    #         merged_assembly=os.path.join(OUTPUTDIR, "Assembly",  "intermediary","h1h2_with_h2h1/merged_with_bridge_singletons_99id.fasta")
    #     output:
    #         final = 'Assembly/%s.assembly.merged.fa' % ASS
    #     log:
    #         os.path.join(OUTPUTDIR, "logs",  "assembly_clean_up_assembly_merge.log")
    #     benchmark:
    #         os.path.join(OUTPUTDIR, "logs",  "assembly_clean_up_assembly_merge_benchmark.tsv")
    #     resources:
    #         runtime = SHORT_RUNTIME,
    #         mem_mb = MAIN_MEMCORE_MB
    #     params:
    #         sample = SAMPLE,
    #         s_dir=os.path.join(OUTPUTDIR,"Assembly","intermediary")
    #     threads:
    #         1
    #     conda:
    #         os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
    #     message:
    #         "Compressing sample {params.sample} intermediary files and creating final merged file."
    #     shell:
    #         """
    #         shopt -s nullglob
    #
    #         # for i in {params.s_dir}/*_formatted.fa {params.s_dir}/join*/*.{{fasta,fasta.clstr}}; do
    #         #   if [[ ${{i}} != {input.merged_assembly} ]]; then
    #         #     echo "Compressing ${{i}}." >> {log} 2>&1
    #         #     gzip -c "${{i}}" > "${{i}}.gz" && rm -f "${{i}}"
    #         #   fi
    #         # done \
    #         # && \
    #         mv {input.merged_assembly} {output.final}
    #         """
