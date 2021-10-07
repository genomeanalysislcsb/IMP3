include:
    "../Assembly/common.smk"

include:
    "../Assembly/mapping.smk"

include:
    "../Assembly/assembly.merge.smk"

if has_sr("mg") and has_sr("mt"):
    include:
        "../Assembly/assemblers/megahit.smk" # TODO: currently the only option for MG/MT assembly
else:
    include:
        "../Assembly/assemblers/%s.smk" % get_assembler()

rule ASSEMBLY_READY:
    input:
        # assembly
        expand(
            "Assembly/{asm}.assembly.merged.{ext}",
            asm=get_omics_assembly(), ext=["fa", "fa.fai", "fa.bed3"]
        ),
        # assembly mapping
        expand(
            "Assembly/{omic}.reads.sorted.{ext}",
            omic=get_omics(), ext=["bam", "bam.bai"]
        ),
        # contig stat.s
        expand(
            "Stats/{omic}/{asm}.assembly.{ext}",
            omic=get_omics(), asm=get_omics_assembly(), ext=["contig_coverage.txt", "contig_depth.txt", "contig_flagstat.txt"]
        ) + expand(
            "Stats/{asm}/{asm}.assembly.{ext}",
            asm=get_omics_assembly(), ext=["length.txt", "gc_content.txt"]
        ),
    output:
        touch("status/assembly.ready")

rule ASSEMBLY:
    input:
        "status/assembly.ready", # not needed (required for the archive anyway) but listed for readability
        "Assembly/intermediary.tar.gz"
    output:
        touch("status/assembly.done")
