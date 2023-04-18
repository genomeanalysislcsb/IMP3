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


# master command
rule ASSEMBLY:
    input:
        # 'Assembly/intermediary.tar.gz',
        'Assembly/mt.assembly.merged.fa',
        'Assembly/mt.assembly.merged.fa.fai',
        'Assembly/mt.assembly.merged.fa.bed3',
        'Assembly/mt.reads.sorted.bam',
        'Assembly/mt.reads.sorted.bam.bai',
        'Stats/mt/mt.assembly.contig_coverage.txt',
        'Stats/mt/mt.assembly.contig_depth.txt',
        'Stats/mt/mt.assembly.contig_flagstat.txt',
        'Stats/mt/mt.assembly.length.txt',
        'Stats/mt/mt.assembly.gc_content.txt'
    output:
        touch('status/assembly.done')
