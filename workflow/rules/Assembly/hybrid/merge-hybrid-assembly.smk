if config["assembly"]["merge"] is None or config["assembly"]["merge"] == "":
    localrules: skip_merge_hybrid_assembly
    rule skip_merge_hybrid_assembly:
        input:
            'Assembly/intermediary/mgmt.%s_hybrid.1.fa' % IMP_ASSEMBLER,
        output:
            'Assembly/mgmt.assembly.merged.fa'
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
            'Assembly/intermediary/mgmt.%s_hybrid.1.fa' % IMP_ASSEMBLER,
            'Assembly/intermediary/mgmt.%s_hybrid.2.fa' % IMP_ASSEMBLER,
        output:
            'Assembly/mgmt.assembly.merged.fa'
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
            i1 = "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % IMP_ASSEMBLER,
            i2 = "Assembly/intermediary/mgmt.%s_hybrid.2.fa" % IMP_ASSEMBLER
        output:
            if1= "Assembly/intermediary/mgmt.%s_hybrid.1_formatted.fa" % IMP_ASSEMBLER,
            if2= "Assembly/intermediary/mgmt.%s_hybrid.2_formatted.fa" % IMP_ASSEMBLER
        log:
            os.path.join(OUTPUTDIR,"logs","assembly_merge_hybrid_format_fastas.log")
        benchmark:
            os.path.join(OUTPUTDIR,"logs","assembly_merge_hybrid_format_fastas_benchmark.tsv")
        resources:
            runtime = "02:00:00",
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

rule quickmerge_join_h1_with_h2:
    input:
        os.path.join(ENVDIR, "qm.ready"),
        assembly_1 = "Assembly/intermediary/mgmt.{assembler}_hybrid.1_formatted.fa",
        assembly_2 = "Assembly/intermediary/mgmt.{assembler}_hybrid.2_formatted.fa"
    output:
        out_dir = directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "{join_prefix}")),
        merged_assembly = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "{join_prefix}/merged_{join_prefix}_rn.fasta")
    log:
        os.path.join(OUTPUTDIR, "logs",  "assembly_{join_prefix}_log")
    benchmark:
        os.path.join(OUTPUTDIR, "logs",  "assembly_{join_prefix}_benchmark.tsv")
    resources:
        runtime = "12:00:00",
        mem_mb = MAIN_MEMCORE_MB
    params:
        assembler=IMP_ASSEMBLER,  # This is the default assembler to be used, can be overridden
        sample = SAMPLE,
        work_dir = OUTPUTDIR,
        min_overlap=config["assembly"]["quickmerge"]["min_overlap"],
        len_cutoff=config["assembly"]["quickmerge"]["len_cutoff"],
        hco=config["assembly"]["quickmerge"]["hco"],
        c=config["assembly"]["quickmerge"]["c"],
        join_prefix="join_h1_with_h2",  # This is the default prefix, can be overridden
        merge_prefix="h1_with_h2"  # This is the default merge prefix, can be overridden
    threads:
        getThreads(128)
    conda:
        os.path.join(ENVDIR, "IMP_assembly_merge.yaml")
    message:
        "Joining sample {params.sample} mgmt hybrid 1 with hybrid 2."
    shell:
        """
        echo -e "File with contigs to join: {input.assembly_1}.\nFile with contigs to use as bridges: {input.assembly_2}." >> {log} 2>&1
        mkdir -p {output.out_dir}

        # We have to change to the output directory because qickmerge will write the output files to the current directory ...
        cd {output.out_dir}

        # Run quickmerge
        time merge_wrapper.py --threads {threads} \
                    --version4 \
                    --prefix "{params.merge_prefix}" \
                    --merging_length_cutoff {params.min_overlap} \
                    --length_cutoff {params.len_cutoff} \
                    --hco {params.hco} \
                    --c {params.c} \
                    {params.work_dir}/{input.assembly_1} \
                    {params.work_dir}/{input.assembly_2} >> {log} 2>&1 \
        
        # Archive the unused quickmerge output files
        parallel --jobs {threads} gzip ::: aln_summary_*.tsv anchor_summary_*.txt param_summary_*.txt {params.merge_prefix}* \
        
        # Rename and format the merged assembly file
        seqkit replace -p .+ \
                        -r "S{params.sample}_{params.merge_prefix}_C{{nr}}" \
                        --nr-width 10 \
                        -o  {output.merged_assembly} \
                        {output.out_dir}/merged_{params.merge_prefix}.fasta \
        
        # Remove the raw merged assembly file
        rm -f {output.out_dir}/merged_{params.merge_prefix}.fasta
        """

use rule quickmerge_join_h1_with_h2 as quickmerge_join_h2_with_h1 with:
    input:
        assembly_1 = "Assembly/intermediary/mgmt.{assembler}_hybrid.2_formatted.fa",
        assembly_2 = "Assembly/intermediary/mgmt.{assembler}_hybrid.1_formatted.fa"
    params:
        join_prefix="join_h2_with_h1",
        merge_prefix="h2_with_h1"
    output:
        out_dir = directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "{join_prefix}")),
        merged_assembly = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "{join_prefix}/merged_{join_prefix}_rn.fasta")
    log:
        os.path.join(OUTPUTDIR, "logs",  "assembly_{join_prefix}_log")
    benchmark:
        os.path.join(OUTPUTDIR, "logs",  "assembly_{join_prefix}_benchmark.tsv")


    rule quickmerge_join_h1h2_with_h2h1:
        input:
            os.path.join(ENVDIR, "qm.ready"),
            assembly_h1h2 = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h1_with_h2/merged_h1_with_h2_rn.fasta"),
            assembly_h2h1 = os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h2_with_h1/merged_h2_with_h1_rn.fasta")
        output:
            out_dir = directory(os.path.join(OUTPUTDIR, "Assembly",  "intermediary", "join_h1h2_with_h2h1")),
            merged_assembly = "Assembly/%s.assembly.merged.fa" %ASS
        log:
            os.path.join(OUTPUTDIR, "logs",  "assembly_join_h1h2_with_h2h1.log")
        benchmark:
            os.path.join(OUTPUTDIR, "logs",  "assembly_join_h1h2_with_h2h1_benchmark.tsv")
        resources:
            runtime = "24:00:00",
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
            # && \
            # rm -rf {output.out_dir}/merged_with_bridge_singletons.fasta {output.out_dir}/merged_h1h2_with_h2h1.fasta \
            #        {output.out_dir}/bridges_h2h1_rn.fasta {output.out_dir}/merged_h1h2_with_h2h1_rn.fasta
            """
