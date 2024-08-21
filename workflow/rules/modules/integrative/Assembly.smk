# include rules for the Assembly part of the workflow

# include needed rules
include:
    '../../Assembly/common/bwa.index.smk'

include:
    '../../Assembly/common/megahit.smk'
# Dont include this, as were already using the hybrid rules
# include:
#     '../../Assembly/common/merge-assembly.smk'

include:
    '../../Assembly/hybrid/%s.hybrid.smk' % IMP_ASSEMBLER

include:
    '../../Assembly/hybrid/extract-unmapped.smk'

include:
    '../../Assembly/hybrid/merge-hybrid-assembly.smk'

include:
    '../../Assembly/common/contig.smk'

include:
    '../../Assembly/common/mapping.smk'

include:
    '../../Assembly/common/contig-length.smk'
include:
    '../../Assembly/common/contig-depth.smk'


# master command
rule ASSEMBLY:
    input:
        'Assembly/%s.assembly.merged.fa' %ASS,
        'Assembly/%s.assembly.merged.fa.fai' %ASS,
        'Assembly/%s.assembly.merged.fa.bed3' %ASS,
        expand('Assembly/{type}.reads.sorted.bam',type=TYPES),
        expand('Assembly/{type}.reads.sorted.bam.bai',type=TYPES),
        expand("Stats/{type}/{ass}.assembly.thresholds.bed.gz",type=TYPES,ass=ASS),
        expand("Stats/{type}/{ass}.assembly.contig_depth.txt",type=TYPES,ass=ASS),
        expand("Stats/{type}/{ass}.assembly.contig_flagstat.txt",type=TYPES,ass=ASS),
        "Stats/%s/%s.assembly.length.txt" % (ASS,ASS),
        "Stats/%s/%s.assembly.gc_content.txt" % (ASS,ASS)
    output:
        touch('status/assembly.done')
