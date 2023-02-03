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
            runtime = "48:00:00",
            mem_mb = BIGMEMCORE_MB
        params:
            sample = SAMPLE
        threads: 1
        conda: os.path.join(ENVDIR, "IMP_assembly.yaml")
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
            runtime = "12:00:00",
            mem_mb = BIGMEMCORE_MB
        params:
            sample = SAMPLE,
            iters = 0,
            minoverlap = 1000
        threads: getThreads(BIGCORENO, bigmem=True)
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        log:
            "logs/assembly_merge_assembly_flye.log"
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
