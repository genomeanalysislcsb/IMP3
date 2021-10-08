# no merging
if not config["assembly"]["merge"]:
    localrules: assembly_merge_skip, assembly_merge_hybrid_skip
    # normal
    rule assembly_merge_skip:
        input:
            "Assembly/intermediary/{asm}.%s_preprocessed.1.fa" % get_assembler() # TODO: rm assembler from the name or add as wildcard
        output:
            temp("Assembly/{asm}.assembly.merged.fa.tmp")
        wildcard_constraints:
            asm=get_omics_assembly(),
        params:
            sample=config["sample"],
        threads: 1
        message:
            "Merge assemblies: skip: use {input} only"
        shell:
            "cp {input} {output}"
    # hybrid
    use rule assembly_merge_skip as assembly_merge_hybrid_skip with:
        input:
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler() # TODO: rm assembler from the name or add as wildcard
        output:
            temp("Assembly/mgmt.assembly.merged.fa.tmp")
        message:
            "Merge hybrid assemblies: skip: use {input} only"

# concat FASTAs
elif config["assembly"]["merge"] == "concat":
    localrules: assembly_merge_concat, assembly_merge_hybrid_concat
    # normal
    rule assembly_merge_concat:
        input:
            "Assembly/intermediary/{asm}.%s_preprocessed.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/{asm}.%s_unmapped.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            temp("Assembly/{asm}.assembly.merged.fa.tmp")
        wildcard_constraints:
            asm=get_omics_assembly(),
        params:
            sample=config["sample"],
        threads: 1
        message:
            "Merge assemblies w/ cat: {input}"
        shell:
            "cat {input} > {output}"
    # hybrid
    use rule assembly_merge_concat as assembly_merge_hybrid_concat with:
        input:
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/mgmt.%s_hybrid.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            temp("Assembly/mgmt.assembly.merged.fa.tmp")
        message:
            "Merge hybrid assemblies w/ cat: {input}"

# use CAP3 for merging
elif config["assembly"]["merge"] == "cap3":
    # normal
    rule assembly_merge_cap3:
        input:
            "Assembly/intermediary/{asm}.%s_preprocessed.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/{asm}.%s_unmapped.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            fa=temp("Assembly/{asm}.assembly.merged.fa.tmp"),
            fa_cat=temp("Assembly/intermediary/{asm}.assembly.cat.fa"),
            fa_contigs=temp("Assembly/intermediary/{asm}.assembly.cat.fa.cap.contigs"),
            fa_singlets=temp("Assembly/intermediary/{asm}.assembly.cat.facap.singlets"),
        log:
            "logs/assembly_merge_cap3.{asm}.log"
        benchmark:
            "benchmark/assembly_merge_cap3.{asm}.tsv"
        wildcard_constraints:
            asm=get_omics_assembly(),
        params:
            sample=config["sample"],
            ident=config["assembly"]["cap3"]["identity"],
            overl=config["assembly"]["cap3"]["overlap"],
        threads: 1
        resources:
            runtime="72:00:00",
            mem=get_mem_big,
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        message:
            "Merge assemblies w/ CAP3: {input}"
        shell:
            # concat input
            "cat {input} > {output.fa_cat} && "
            # run CAP3
            "cap3 {output.fa_cat} -p {params.ident} -o {params.overl} &> {log} && "
            # concat assembled contigs and singletons | rename the contigs
            "cat {output.fa_contigs} {output.fa_singlets} | awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > {output.fa}"
    
    # hybrid
    use rule assembly_merge_cap3 as assembly_merge_hybrid_cap3 with:
        input:
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/mgmt.%s_hybrid.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            fa=temp("Assembly/mgmt.assembly.merged.fa.tmp"),
            fa_cat=temp("Assembly/intermediary/mgmt.assembly.cat.fa"),
            fa_contigs=temp("Assembly/intermediary/mgmt.assembly.cat.fa.cap.contigs"),
            fa_singlets=temp("Assembly/intermediary/mgmt.assembly.cat.facap.singlets"),
        log:
            "logs/assembly_merge_hybrid_cap3.log"
        benchmark:
            "benchmark/assembly_merge_hybrid_cap3.tsv"
        message:
            "Merge hybrid assemblies w/ CAP3: {input}"

elif config["assembly"]["merge"] == "flye":
    # normal
    rule assembly_merge_flye:
        input:
            i1="Assembly/intermediary/{asm}.%s_preprocessed.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            i2="Assembly/intermediary/{asm}.%s_unmapped.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
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
            final=temp("Assembly/{asm}.assembly.merged.fa.tmp"),
        log:
            "logs/assembly_merge_flye.{asm}.log"
        benchmark:
            "benchmark/assembly_merge_flye.{asm}.tsv"
        wildcard_constraints:
            asm=get_omics_assembly(),
        threads: get_cores_big()
        resources:
            runtime="24:00:00",
            mem=get_mem_big,
        params:
            iters=0,
            minoverlap=1000
        conda:
            os.path.join(ENVDIR, "IMP_assembly.yaml")
        message:
            "Merge assemblies w/ Flye: {input}"
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
            cat {output.flye} {output.ufna} > {output.final}
            ) &> {log}
            """
    # hybrid
    use rule assembly_merge_flye as assembly_merge_hybrid_flye with:
        input:
            i1="Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            i2="Assembly/intermediary/mgmt.%s_hybrid.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
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
            final=temp("Assembly/mgmt.assembly.merged.fa.tmp"),
        log:
            "logs/assembly_merge_hybrid_flye.log"
        benchmark:
            "benchmark/assembly_merge_hybrid_flye.tsv"
        message:
            "Merge hybrid assemblies w/ Flye: {input}"

localrules: rename_contigs
rule assembly_rename_contigs:
    input:
        "Assembly/{asm}.assembly.merged.fa.tmp"
    output:
        "Assembly/{asm}.assembly.merged.fa"
    # wildcard_constraints:
        # TODO
    params:
        sample=config["sample"]
    message:
        "Rename contigs: {input}"
    shell:
        "awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {input} > {output}"