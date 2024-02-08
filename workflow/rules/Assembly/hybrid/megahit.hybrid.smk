rule megahit_hybrid_assembly_1:
    input:
        mg_r1 = 'Preprocessing/mg.r1.preprocessed.fq',
        mg_r2 = 'Preprocessing/mg.r2.preprocessed.fq',
        mg_se = 'Preprocessing/mg.se.preprocessed.fq',
        mt_r1 = 'Preprocessing/mt.r1.preprocessed.fq',
        mt_r2 = 'Preprocessing/mt.r2.preprocessed.fq',
        mt_se = 'Preprocessing/mt.se.preprocessed.fq',
        mt_ass_pre = 'Assembly/intermediary/mt.megahit_preprocessed.1/final.contigs.fa',
        mt_unmapped = 'Assembly/intermediary/mt.megahit_unmapped.2/final.contigs.fa',
    output:
        contigs_final = 'Assembly/intermediary/mgmt.megahit_hybrid.1/final.contigs.fa',
        contigs_hyb = 'Assembly/intermediary/mgmt.megahit_hybrid.1.fa',
        out_dir = directory('Assembly/intermediary/mgmt.megahit_hybrid.1')
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        max_mem = int(MAIN_BIGMEMTOTAL * byte_dict[config['cluster_specs']["mem_unit"]]),
        kmer_list = config["assembly"]["kmer_list"],
        # Wrap se inputs as a parameter for use in the shell command.
        single = lambda wildcards, input: ','.join([input.mg_se, input.mt_se, input.mt_ass_pre, input.mt_unmapped]),
    conda:
        ENVDIR + "/IMP_assembly.yaml"
    log:
        "logs/assembly/megahit_hybrid_assembly_1.log"
    message:
        "megahit_hybrid_assembly_1: Performing hybrid assembly 1 from preprocessed reads using MEGAHIT"
    shell:
        """
        if [ -d "{output.out_dir}" ]; then
            rm -rf {output.out_dir}
        fi

        megahit -1 {input.mg_r1},{input.mt_r1} \
                -2 {input.mg_r2},{input.mt_r2} \
                -r {params.single} \
                --out-dir {output.out_dir} \
                --k-list {params.kmer_list} \
                --num-cpu-threads {threads} \
                --memory {params.max_mem} \
                --mem-flag 1 > {log} 2>&1
        ln -fs $(echo {output.contigs_final} | cut -f 3,4 -d /) {output.contigs_hyb} && touch -h {output.contigs_hyb}
        """

use rule megahit_hybrid_assembly_1 as megahit_hybrid_assembly_2 with:
    input:
        mg_r1 = 'Assembly/intermediary/mgmt.r1.megahit_hybrid.mt.unmapped.fq',
        mg_r2 = 'Assembly/intermediary/mgmt.r2.megahit_hybrid.mt.unmapped.fq',
        mg_se = 'Assembly/intermediary/mgmt.se.megahit_hybrid.mt.unmapped.fq',
        mt_r1 = 'Assembly/intermediary/mgmt.r1.megahit_hybrid.mg.unmapped.fq',
        mt_r2 = 'Assembly/intermediary/mgmt.r2.megahit_hybrid.mg.unmapped.fq',
        mt_se = 'Assembly/intermediary/mgmt.se.megahit_hybrid.mg.unmapped.fq',
    output:
        contigs_final = 'Assembly/intermediary/mgmt.megahit_hybrid.2/final.contigs.fa',
        contigs_hyb = 'Assembly/intermediary/mgmt.megahit_hybrid.2.fa',
        out_dir = directory('Assembly/intermediary/mgmt.megahit_hybrid.2')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        max_mem = int(MAIN_BIGMEMTOTAL * byte_dict[config['cluster_specs']["mem_unit"]]),
        kmer_list = config["assembly"]["kmer_list"],
        single = lambda wildcards, input: ','.join([input.mg_se, input.mt_se]),
    conda:
        ENVDIR + "/IMP_assembly.yaml"
    log:
        "logs/assembly/megahit_hybrid_assembly_2.log"
    message:
        "megahit_hybrid_assembly_2: Performing hybrid assembly 2 from preprocessed reads using MEGAHIT"
