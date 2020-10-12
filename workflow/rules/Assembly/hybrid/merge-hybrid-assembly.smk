if config["assembly"]["merge"]:
    rule merge_hybrid_assembly_cap3:
        input:
            'Assembly/intermediary/mgmt.%s_hybrid.1.fa' % IMP_ASSEMBLER,
            'Assembly/intermediary/mgmt.%s_hybrid.2.fa' % IMP_ASSEMBLER,
        output:
            'Assembly/mgmt.assembly.merged.fa'
        params:
            runtime = "72:00:00",
            mem = BIGMEMCORE,
            sample = SAMPLE
        threads: 1
        conda: ENVDIR + "/IMP_assembly.yaml"
        log: "assembly_merge_hybrid_assembly_cap3.log"
        message: "merge_hybrid_assembly_cap3: Merging hybrid assemblies."
        shell:
            """
            NAME_fin=Assembly/mgmt.assembly
            NAME=Assembly/intermediary/mgmt.assembly
            cat {input} > $NAME.cat.fa

            # Options should be after input file!
            cap3 $NAME.cat.fa -p {config[assembly][cap3][identity]} -o {config[assembly][cap3][overlap]} > {log} 2>&1

            # Concatenate assembled contigs, singletons and rename the contigs
            cat $NAME.cat.fa.cap.contigs $NAME.cat.fa.cap.singlets | \
            awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > $NAME_fin.merged.fa
            """
else:
    localrules: skip_merge_hybrid_assembly_cap3
    rule skip_merge_hybrid_assembly_cap3:
        input:
            'Assembly/intermediary/mgmt.%s_hybrid.1.fa' % IMP_ASSEMBLER,
            'Assembly/intermediary/mgmt.%s_hybrid.2.fa' % IMP_ASSEMBLER,
        output:
            'Assembly/mgmt.assembly.merged.fa'
        params:
            sample = SAMPLE
        threads: 1
        message:
            "skip_merge_hybrid_assembly_cap3: Not merging hybrid assemblies, copying {input[0]}."
        shell:
            "awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {input[0]} > {output}"
