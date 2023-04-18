# include rules for the Assembly part of the workflow

# include needed rules
include:
    '../../Assembly/common/bwa.index.smk'

include:
    '../../Assembly/common/megahit.smk'
    
include:
    '../../Assembly/common/merge-assembly.smk'

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
        # 'Assembly/intermediary.tar.gz',
        'Assembly/%s.assembly.merged.fa' %ASS,
        'Assembly/%s.assembly.merged.fa.fai' %ASS,
        'Assembly/%s.assembly.merged.fa.bed3' %ASS,
        'Assembly/mt.reads.sorted.bam',
        'Assembly/mt.reads.sorted.bam.bai',
        'Assembly/mg.reads.sorted.bam',
        'Assembly/mg.reads.sorted.bam.bai',
        expand("Stats/{type}/{ass}.assembly.contig_coverage.txt",type=TYPES,ass=ASS),
        expand("Stats/{type}/{ass}.assembly.contig_depth.txt",type=TYPES,ass=ASS),
        expand("Stats/{type}/{ass}.assembly.contig_flagstat.txt",type=TYPES,ass=ASS),
        "Stats/%s/%s.assembly.length.txt" % (ASS,ASS),
        "Stats/%s/%s.assembly.gc_content.txt" % (ASS,ASS)
    output:
        touch('status/assembly.done')
