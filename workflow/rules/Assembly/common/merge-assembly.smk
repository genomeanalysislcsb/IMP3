if config["assembly"]["merge"]:
    rule merge_assembly_cap3:
        input:
            'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            'Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            'Assembly/%s.assembly.merged.fa' % ASS,
        params:
            runtime = "72:00:00",
            mem = BIGMEMCORE,
            sample = SAMPLE
        threads: 1
        conda: ENVDIR + "/IMP_assembly.yaml"
        log: "logs/assembly_merge_assembly_cap3.log"
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
else:
    localrules: skip_merge_assembly_cap3
    rule skip_merge_assembly_cap3:
        input:
            'Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            'Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            'Assembly/%s.assembly.merged.fa' % ASS,
        params:
            sample = SAMPLE
        threads: 1
        log:
            "logs/assembly_skip_merge_assembly_cap3.log"
        message:
            "skip_merge_assembly_cap3: Not merging assemblies, copying {input[0]}."
        shell:
            "awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {input[0]} > {output}"
