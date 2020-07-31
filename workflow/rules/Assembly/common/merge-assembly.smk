if config["assembly"]["merger"] == "cap3":
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
elif config["assembly"]["merger"] == "flye":
    rule merge_assembly_flye:
        input:
            i1='Assembly/intermediary/%s.%s_preprocessed.1.fa' % (ASS, IMP_ASSEMBLER),
            i2='Assembly/intermediary/%s.%s_unmapped.2.fa' % (ASS, IMP_ASSEMBLER),
        output:
            # temp files (processed input)
            i1tmp=temp("Assembly/intermediary/common_merge_flye/input.1.fasta"),
            i2tmp=temp("Assembly/intermediary/common_merge_flye/input.2.fasta"),
            # intermediate files
            flye="Assembly/intermediary/common_merge_flye/assembly.fasta",
            # final files
            final='Assembly/%s.assembly.merged.fa' % ASS,
        params:
            runtime = "12:00:00",
            mem = BIGMEMCORE,
            sample = SAMPLE,
            iters = 0
        threads: BIGCORENO
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        log:
            "logs/assembly_merge_assembly_flye.log"
        message:
            "merge_assembly_flye: Merging assemblies."
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
