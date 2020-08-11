if config["assembly"]["merger"] == "cap3":
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
elif config["assembly"]["merger"] == "flye":
    rule merge_hybrid_assembly_flye:
        input:
            i1="Assembly/intermediary/mgmt.%s_hybrid.1.fa" % IMP_ASSEMBLER,
            i2="Assembly/intermediary/mgmt.%s_hybrid.2.fa" % IMP_ASSEMBLER,
        output:
            # temp files (processed input)
            i1tmp=temp("Assembly/intermediary/hybrid_merge_flye/input.1.fasta"),
            i2tmp=temp("Assembly/intermediary/hybrid_merge_flye/input.2.fasta"),
            # intermediate files
            flye="Assembly/intermediary/hybrid_merge_flye/assembly.fasta",
            # final files
            final="Assembly/mgmt.assembly.merged.fa"
        params:
            runtime = "24:00:00",
            mem = BIGMEMCORE,
            sample = SAMPLE,
            iters = 0
        threads: BIGCORENO
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        log:
            "logs/assembly_merge_hybrid_assembly_flye.log"
        message:
            "merge_hybrid_assembly_flye: Merging hybrid assemblies."
        shell:
            """
            # genome size estimation in Mb
            genomesize=$(grep -v "^>" {input[0]} | awk '{{sum += length}} END {{printf \"%.0f\\n\", sum/1e+6}}')
            # rename contigs in input files
            awk '/^>/{{print \">input1_contig_\" ++i; next}}{{print}}' {input.i1} > {output.i1tmp}
            awk '/^>/{{print \">input2_contig_\" ++i; next}}{{print}}' {input.i2} > {output.i2tmp}
            # consensus assembly
            flye --subassemblies {output.i1tmp} {output.i2tmp} -o $(dirname {output.flye}) -g \"${{genomesize}}m\" -t {threads} -i {params.iters}
            # rename contigs
            awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {output.flye} > {output.final}
            """
