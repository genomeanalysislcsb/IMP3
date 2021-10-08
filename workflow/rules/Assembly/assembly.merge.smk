# no merging
if not config["assembly"]["merge"]:
    localrules: merge_assembly_skip, merge_hybrid_assembly_skip
    # normal
    rule merge_assembly_skip:
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
    use rule merge_assembly_skip as merge_hybrid_assembly_skip with:
        input:
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler() # TODO: rm assembler from the name or add as wildcard
        output:
            temp("Assembly/mgmt.assembly.merged.fa.tmp")
        message:
            "Merge hybrid assemblies: skip: use {input} only"

# concat FASTAs
elif config["assembly"]["merge"] == "concat":
    localrules: merge_assembly_concat, merge_hybrid_assembly_concat
    # normal
    rule merge_assembly_concat:
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
            "Merge assemblies: concat: {input}"
        shell:
            "cat {input} > {output}"
    # hybrid
    use rule merge_assembly_concat as merge_hybrid_assembly_concat with:
        input:
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/mgmt.%s_hybrid.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            temp("Assembly/mgmt.assembly.merged.fa.tmp")
        message:
            "Merge hybrid assemblies: concat: {input}"

# use CAP3 for merging
elif config["assembly"]["merge"] == "cap3":
    # normal
    rule merge_assembly_cap3:
        input:
            "Assembly/intermediary/{asm}.%s_preprocessed.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/{asm}.%s_unmapped.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            fa=temp("Assembly/{asm}.assembly.merged.fa.tmp"),
            fa_cat=temp("Assembly/intermediary/{asm}.assembly.cat.fa"),
            fa_contigs=temp("Assembly/intermediary/{asm}.assembly.cat.fa.cap.contigs"),
            fa_singlets=temp("Assembly/intermediary/{asm}.assembly.cat.facap.singlets"),
        log:
            "logs/assembly_merging.cap3.{asm}.log"
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
            "Merge assemblies: CAP3: {input}"
        shell:
            # concat input
            "cat {input} > {output.fa_cat} && "
            # run CAP3
            "cap3 {output.fa_cat} -p {params.ident} -o {params.overl} &> {log} && "
            # concat assembled contigs and singletons | rename the contigs
            "cat {output.fa_contigs} {output.fa_singlets} | awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' > {output.fa}"
    
    # hybrid
    use rule merge_assembly_cap3 as merge_hybrid_assembly_cap3 with:
        input:
            "Assembly/intermediary/mgmt.%s_hybrid.1.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
            "Assembly/intermediary/mgmt.%s_hybrid.2.fa" % get_assembler(), # TODO: rm assembler from the name or add as wildcard
        output:
            fa=temp("Assembly/mgmt.assembly.merged.fa.tmp"),
            fa_cat=temp("Assembly/intermediary/mgmt.assembly.cat.fa"),
            fa_contigs=temp("Assembly/intermediary/mgmt.assembly.cat.fa.cap.contigs"),
            fa_singlets=temp("Assembly/intermediary/mgmt.assembly.cat.facap.singlets"),
        log:
            "logs/assembly_merging.hybrid.cap3.log"
        message:
            "Merge hybrid assemblies: CAP3: {input}"

elif config["assembly"]["merge"] == "flye":
    # normal
    rule merge_assembly_flye:
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
            "logs/assembly_merge_assembly_flye.log"
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
            "Merge assemblies: Flye: {input}"
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
    use rule merge_assembly_flye as merge_hybrid_assembly_flye with:
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
            "logs/assembly_merge_hybrid_assembly_flye.log"
        message:
            "merge_hybrid_assembly_flye: Merging hybrid assemblies."

localrules: rename_contigs
rule rename_contigs:
    input:
        "Assembly/{asm}.assembly.merged.fa.tmp"
    output:
        "Assembly/{asm}.assembly.merged.fa"
    # wildcard_constraints:
        # TODO
    params:
        sample=config["sample"]
    shell:
        "awk '/^>/{{print \">{params.sample}_contig_\" ++i; next}}{{print}}' {input} > {output}"