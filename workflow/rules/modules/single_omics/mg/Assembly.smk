# include rules for the Assembly part of the workflow

# include needed rules
include:
    '../../../Assembly/common/bwa.index.smk'

include:
    '../../../Assembly/common/%s.smk' % IMP_ASSEMBLER

include:
    '../../../Assembly/common/merge-assembly.smk'

include:
    '../../../Assembly/single-omic/extract-unmapped.smk'

include:
    '../../../Assembly/common/contig.smk'

include:
    '../../../Assembly/common/mapping.smk'

include:
    '../../../Assembly/common/contig-length.smk'

include:
    '../../../Assembly/common/contig-depth.smk'


master_ass_input = [
        # 'Assembly/intermediary.tar.gz',
        'Assembly/mg.assembly.merged.fa',
        'Assembly/mg.assembly.merged.fa.fai',
        'Assembly/mg.assembly.merged.fa.bed3',
        expand('Assembly/{type}.reads.sorted.bam',type=TYPES),
        expand('Assembly/{type}.reads.sorted.bam.bai',type=TYPES),
        expand("Stats/{type}/{ass}.assembly.contig_coverage.txt",type=TYPES,ass=ASS),
        expand("Stats/{type}/{ass}.assembly.contig_depth.txt",type=TYPES,ass=ASS),
        expand("Stats/{type}/{ass}.assembly.contig_flagstat.txt",type=TYPES,ass=ASS),
        "Stats/%s/%s.assembly.length.txt" % (ASS,ASS),
        "Stats/%s/%s.assembly.gc_content.txt" % (ASS,ASS)]

# master command
rule ASSEMBLY:
    input:
        master_ass_input
    output:
        touch('status/assembly.done')
