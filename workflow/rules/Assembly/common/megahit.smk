MEGAHIT_ASSEMBLY_SHELL = """
                         if [ -d "{output.out_dir}" ]; then
                           rm -rf {output.out_dir}
                         fi

                         megahit -1 {input.r1} \
                                 -2 {input.r2} \
                                 -r {input.se} \
                                 --out-dir {output.out_dir} \
                                 --k-list {params.kmer_list} \
                                 --num-cpu-threads {threads} \
                                 -m {params.max_mem} \
                                 --mem-flag 1 > {log} 2>&1
                         ln -fs $(echo {output.contigs_final} | cut -f 3,4 -d /) {output.contigs} && touch -h {output.contigs}
                         """

rule megahit_assembly_from_preprocessing:
    input:
        r1 = 'Preprocessing/{type}.r1.preprocessed.fq',
        r2 = 'Preprocessing/{type}.r2.preprocessed.fq',
        se = 'Preprocessing/{type}.se.preprocessed.fq'
    output:
        out_dir = directory('Assembly/intermediary/{type}.megahit_preprocessed.1'),
        contigs_final = 'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa',
        contigs = 'Assembly/intermediary/{type}.megahit_preprocessed.1.fa'
    resources:
        runtime = "24:00:00",
        mem_mb = BIGMEMCORE_MB
    threads:
        getThreads(BIGCORENO, bigmem=True)
    params:
        max_mem = int(BIGMEMTOTAL * byte_dict[config["mem"]["mem_unit"]]),
        kmer_list = config["assembly"]["kmer_list"]
    conda:
        ENVDIR + "/IMP_assembly.yaml"
    log:
        "logs/assembly_megahit_assembly_from_preprocessing.{type}.log"
    message:
        "megahit_assembly_from_preprocessing: Performing {wildcards.type} assembly step 1 from preprocessed reads using MEGAHIT"
    shell:
        MEGAHIT_ASSEMBLY_SHELL

rule megahit_assembly_from_unmapped:
    input:
        r1 = 'Assembly/intermediary/{type}.r1.unmapped.fq',
        r2 = 'Assembly/intermediary/{type}.r2.unmapped.fq',
        se = 'Assembly/intermediary/{type}.se.unmapped.fq'
    output:
        out_dir = directory('Assembly/intermediary/{type}.megahit_unmapped.2'),
        contigs_final = 'Assembly/intermediary/{type}.megahit_unmapped.2/final.contigs.fa',
        contigs = 'Assembly/intermediary/{type}.megahit_unmapped.2.fa'
    resources:
        runtime = "24:00:00",
        mem_mb = BIGMEMCORE_MB
    threads:
        getThreads(BIGCORENO, bigmem=True)
    params:
        max_mem = int(BIGMEMTOTAL * byte_dict[config["mem"]["mem_unit"]]),
        kmer_list = config["assembly"]["kmer_list"]
    conda:
        ENVDIR + "/IMP_assembly.yaml"
    log:
        "logs/assembly_megahit_assembly_from_unmapped.{type}.log"
    message:
        "megahit_assembly_from_unmapped: Performing {wildcards.type} assembly step 2 from unmapped reads using MEGAHIT"
    shell:
        MEGAHIT_ASSEMBLY_SHELL
