# TODO

include:
    '../Assembly/common.smk'

include:
    '../Assembly/mapping.smk'

include:
    '../Assembly/assembly.merge.smk'

if has_sr("mg") and has_sr("mt"):
    include:
        '../Assembly/assemblers/megahit.smk' # TODO: currently the only option for MG/MT assembly

    rule ASSEMBLY:
        input:
            # 'Assembly/intermediary.tar.gz',
            'Assembly/%s.assembly.merged.fa' % get_omics_assembly(),
            # 'Assembly/%s.assembly.merged.fa.fai' %ASS,
            # 'Assembly/%s.assembly.merged.fa.bed3' %ASS,
            # 'Assembly/mt.reads.sorted.bam',
            # 'Assembly/mt.reads.sorted.bam.bai',
            # 'Assembly/mg.reads.sorted.bam',
            # 'Assembly/mg.reads.sorted.bam.bai',
            # expand("Stats/{type}/{ass}.assembly.contig_coverage.txt",type=TYPES,ass=ASS),
            # expand("Stats/{type}/{ass}.assembly.contig_depth.txt",type=TYPES,ass=ASS),
            # expand("Stats/{type}/{ass}.assembly.contig_flagstat.txt",type=TYPES,ass=ASS),
            # "Stats/%s/%s.assembly.length.txt" % (ASS,ASS),
            # "Stats/%s/%s.assembly.gc_content.txt" % (ASS,ASS)
        output:
            touch('status/assembly.done')
# MG only
elif has_sr("mg"):
    include:
        '../Assembly/assemblers/%s.smk' % get_assembler()

    rule ASSEMBLY:
        input:
            # 'Assembly/intermediary.tar.gz',
            'Assembly/mg.assembly.merged.fa',
            # 'Assembly/mg.assembly.merged.fa.fai',
            # 'Assembly/mg.assembly.merged.fa.bed3',
            # expand('Assembly/{type}.reads.sorted.bam',type=TYPES),
            # expand('Assembly/{type}.reads.sorted.bam.bai',type=TYPES),
            # expand("Stats/{type}/{ass}.assembly.contig_coverage.txt",type=TYPES,ass=ASS),
            # expand("Stats/{type}/{ass}.assembly.contig_depth.txt",type=TYPES,ass=ASS),
            # expand("Stats/{type}/{ass}.assembly.contig_flagstat.txt",type=TYPES,ass=ASS),
            # "Stats/%s/%s.assembly.length.txt" % (ASS,ASS),
            # "Stats/%s/%s.assembly.gc_content.txt" % (ASS,ASS)
        output:
            touch('status/assembly.done')
# MT only
elif has_sr("mt"):
    include:
        '../Assembly/assemblers/%s.smk' % get_assembler()
    
    rule ASSEMBLY:
        input:
            # 'Assembly/intermediary.tar.gz',
            'Assembly/mt.assembly.merged.fa',
            # 'Assembly/mt.assembly.merged.fa.fai',
            # 'Assembly/mt.assembly.merged.fa.bed3',
            # 'Assembly/mt.reads.sorted.bam',
            # 'Assembly/mt.reads.sorted.bam.bai',
            # 'Stats/mt/mt.assembly.contig_coverage.txt',
            # 'Stats/mt/mt.assembly.contig_depth.txt',
            # 'Stats/mt/mt.assembly.contig_flagstat.txt',
            # 'Stats/mt/mt.assembly.length.txt',
            # 'Stats/mt/mt.assembly.gc_content.txt'
        output:
            touch('status/assembly.done')
