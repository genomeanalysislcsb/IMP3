elif config["assembly"]["merge"] is None or config["assembly"]["merge"] == "":
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
elif config["assembly"]["merge"] == "cap3":
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
elif config["assembly"]["merge"] == "flye":
    rule merge_hybrid_assembly_flye:
        input:
            i1="Assembly/intermediary/mgmt.%s_hybrid.1.fa" % IMP_ASSEMBLER,
            i2="Assembly/intermediary/mgmt.%s_hybrid.2.fa" % IMP_ASSEMBLER,
        output:
            # temp files (processed input)
            i1=temp("Assembly/intermediary/hybrid_merge_flye/input.1.fasta"),
            i2=temp("Assembly/intermediary/hybrid_merge_flye/input.2.fasta"),
            i12=temp("Assembly/intermediary/hybrid_merge_flye/input.fasta"),
            uids=temp("Assembly/intermediary/hybrid_merge_flye/unmapped.ids"),
            ufna=temp("Assembly/intermediary/hybrid_merge_flye/unmapped.fasta"),
            # intermediate files
            flye="Assembly/intermediary/hybrid_merge_flye/assembly.fasta",
            # final files
            final="Assembly/mgmt.assembly.merged.fa"
        params:
            runtime = "24:00:00",
            mem = BIGMEMCORE,
            sample = SAMPLE,
            iters = 0,
            minoverlap = 1000,
        threads: BIGCORENO
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        log:
            "logs/assembly_merge_hybrid_assembly_flye.log"
        message:
            "merge_hybrid_assembly_flye: Merging hybrid assemblies."
        shell:
            """
            (
            # rename contigs in input files
            awk '/^>/{{print \">input1_contig_\" ++i; next}}{{print}}' {input.i1} > {output.i1}
            awk '/^>/{{print \">input2_contig_\" ++i; next}}{{print}}' {input.i2} > {output.i2}
            cat {output.i1} {output.i2} > {output.i12}
            # consensus assembly
            flye --subassemblies {output.i1tmp} {output.i2tmp} -o $(dirname {output.flye}) -t {threads} -i {params.iters} --min-overlap {params.minoverlap} 1>&2
            # align input contigs to to assembly | get unmapped contigs | get contig IDs
            flye-minimap2 -ax asm5 {output.flye} {output.i12} | samtools view -f 4 - | cut -f1 > {output.uids}
            # filter input contigs by contig ID
            seqkit grep -w 0 --by-name -f {output.uids} {output.i12} -o {output.ufna}
            # add unmapped contigs | rename contigs
            cat {output.flye} {output.ufna} | awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > {output.final}
            ) &> {log}
            """
