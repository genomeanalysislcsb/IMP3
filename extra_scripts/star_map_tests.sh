
micromamba activate /mnt/isilon/projects/ecosystem_biology/local_tools/envs/imp3/ddc58aa5e3d04a1cbcdc87e209214517_

star_bin="/mnt/isilon/projects/ecosystem_biology/local_tools/envs/imp3/ddc58aa5e3d04a1cbcdc87e209214517_/bin/STAR-avx2"

ulimit -n 10000
# mt.reads.paired.sorted

${star_bin} --runThreadN 32 \
            --runMode alignReads \
            --genomeDir Assembly_test \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 150 \
            --alignIntronMax 1 \
            --outFileNamePrefix Assembly_test/ \
            --outSAMtype BAM \
                         SortedByCoordinate \
            --outSAMheaderHD "@RG\tID:imp3_test\tSM:mt\tPL:platform\tLB:library" \
            --outBAMcompression -1 \
            --outReadsUnmapped Fastx \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --readFilesIn Preprocessing/mt.r1.preprocessed.fq \
                          Preprocessing/mt.r2.preprocessed.fq