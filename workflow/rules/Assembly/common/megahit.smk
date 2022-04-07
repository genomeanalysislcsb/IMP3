MEGAHIT_ASSEMBLY_SHELL = """
if [ -d "{output[0]}" ]; then
    rm -rf {output[0]}
fi

megahit -1 {input[0]} \
 -2 {input[1]} \
 -r {input[2]} \
 -o {output[0]} \
 --k-min {config[assembly][mink]} \
 --k-max {config[assembly][maxk]} \
 --k-step {config[assembly][step]} \
 --bubble-level 0 \
 -t {threads} --cpu-only \
 -m {params.max_mem} \
 --mem-flag 1 > {log} 2>&1
ln -fs $(echo {output[1]} | cut -f 3,4 -d /) {output[2]} && touch -h {output[2]}
"""

rule megahit_assembly_from_preprocessing:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq'
    output:
        directory('Assembly/intermediary/{type}.megahit_preprocessed.1'),
        'Assembly/intermediary/{type}.megahit_preprocessed.1/final.contigs.fa',
        'Assembly/intermediary/{type}.megahit_preprocessed.1.fa'
    resources:
        runtime = "120:00:00",
        mem = BIGMEMCORE
    threads:
        getThreads(BIGCORENO)
    params:
        max_mem = int(BIGMEMTOTAL * byte_dict[config["mem"]["mem_unit"]])
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
        'Assembly/intermediary/{type}.r1.unmapped.fq',
        'Assembly/intermediary/{type}.r2.unmapped.fq',
        'Assembly/intermediary/{type}.se.unmapped.fq'
    output:
        directory('Assembly/intermediary/{type}.megahit_unmapped.2'),
        'Assembly/intermediary/{type}.megahit_unmapped.2/final.contigs.fa',
        'Assembly/intermediary/{type}.megahit_unmapped.2.fa'
    resources:
        runtime = "120:00:00",
        mem = BIGMEMCORE
    threads:
        getThreads(BIGCORENO)
    params:
        max_mem = int(BIGMEMTOTAL * byte_dict[config["mem"]["mem_unit"]])
    conda:
        ENVDIR + "/IMP_assembly.yaml"
    log:
        "logs/assembly_megahit_assembly_from_unmapped.{type}.log"
    message:
        "megahit_assembly_from_unmapped: Performing {wildcards.type} assembly step 2 from unmapped reads using MEGAHIT"
    shell:
        MEGAHIT_ASSEMBLY_SHELL
