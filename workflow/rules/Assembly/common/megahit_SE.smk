rule megahit_assembly_from_preprocessing_se:
    input:
        se = 'Preprocessing/{type}.se.preprocessed.fq'
    output:
        out_dir = directory('Assembly/intermediary/{type}.megahit_preprocessed.1'),
        contigs_final = 'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa',
        contigs = 'Assembly/intermediary/{type}.megahit_preprocessed.1.fa'
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        max_mem = int(MAIN_BIGMEMTOTAL * byte_dict[config['cluster_specs']["mem_unit"]]),
        kmer_list = config["assembly"]["kmer_list"]
    conda:
        ENVDIR + "/IMP_assembly.yaml"
    log:
        "logs/assembly/megahit_assembly_from_preprocessing.{type}.log"
    message:
        "megahit_assembly_from_preprocessing: Performing {wildcards.type} assembly step 1 from preprocessed reads using MEGAHIT"
    shell:
        """
        if [ -d "{output.out_dir}" ]; then
          rm -rf {output.out_dir}
        fi

        megahit -r {input.se} \
                --out-dir {output.out_dir} \
                --k-list {params.kmer_list} \
                --num-cpu-threads {threads} \
                -m {params.max_mem} \
                --mem-flag 1 > {log} 2>&1
        ln -fs $(echo {output.contigs_final} | cut -f 3,4 -d /) {output.contigs} && touch -h {output.contigs}
        """

use rule megahit_assembly_from_preprocessing_se as megahit_assembly_from_unmapped_se with:
    input:
        se = 'Assembly/intermediary/{type}.se.unmapped.fq'
    output:
        out_dir = directory('Assembly/intermediary/{type}.megahit_unmapped.2'),
        contigs_final = 'Assembly/intermediary/{type}.megahit_unmapped.2/final.contigs.fa',
        contigs = 'Assembly/intermediary/{type}.megahit_unmapped.2.fa'
    log:
        "logs/assembly/megahit_assembly_from_unmapped.{type}.log"
    message:
        "megahit_assembly_from_unmapped: Performing {wildcards.type} assembly step 2 from unmapped reads using MEGAHIT"
