METASPADES_ASSEMBLY_SHELL = """
if [ -d "{output[0]}" ]; then
    rm -rf {output[0]}
fi
spades.py --meta \
 --pe1-1 {input[0]} \
 --pe1-2 {input[1]} \
 --pe1-s {input[2]} \
  -t {threads} \
  -m {BIGMEMTOTAL} \
  -k {KMER_STEPS} \
  {LONG_READ_ARG} \
  -o {output[0]} > {log} 2>&1
ln -fs  {output[1]} {output[2]}
"""

kmersteps=range(config['assembly']['mink'],config['assembly']['maxk']+1,config['assembly']['step'])
KMER_STEPS =",".join(map(str,kmersteps)) 

LONG_READ_ARG = ""
if LONG:
    LONG_READ_ARG = "--" + config['raws']['LongReadTech'] + " Assembly/long_reads.fastq"
    rule metaspades_assembly_from_preprocessing:
        input:
            'Preprocessing/{type}.r1.preprocessed.fq',
            'Preprocessing/{type}.r2.preprocessed.fq',
            'Preprocessing/{type}.se.preprocessed.fq',
            'Assembly/long_reads.fastq'
        output:
            directory('Assembly/intermediary/{type}.metaspades_preprocessed.1'),
            'Assembly/intermediary/{type}.metaspades_preprocessed.1/contigs.fa',
            'Assembly/intermediary/{type}.metaspades_preprocessed.1.fa'
        params:
            runtime = "120:00:00",
            mem = BIGMEMCORE
        threads: BIGCORENO
        conda: ROOT_DIR + "/envs/IMP_assembly.yaml"
        log: "logs/assembly_metaspades_assembly_from_preprocessing.{type}.log"
        message: "metaspades_assembly_from_preprocessing: Performing {wildcards.type} assembly step 1 from preprocessed reads using MetaSpades"
        shell:
            METASPADES_ASSEMBLY_SHELL
else:
    rule metaspades_assembly_from_preprocessing:
        input:
            'Preprocessing/{type}.r1.preprocessed.fq',
            'Preprocessing/{type}.r2.preprocessed.fq',
            'Preprocessing/{type}.se.preprocessed.fq'
        output:
            directory('Assembly/intermediary/{type}.metaspades_preprocessed.1'),
            'Assembly/intermediary/{type}.metaspades_preprocessed.1/contigs.fa',
            'Assembly/intermediary/{type}.metaspades_preprocessed.1.fa'
        params:
            runtime = "120:00:00",
            mem = BIGMEMCORE
        threads: BIGCORENO
        conda: ROOT_DIR + "/envs/IMP_assembly.yaml"
        log: "logs/assembly_metaspades_assembly_from_preprocessing.{type}.log"
        message: "metaspades_assembly_from_preprocessing: Performing {wildcards.type} assembly step 1 from preprocessed reads using MetaSpades"
        shell:
            METASPADES_ASSEMBLY_SHELL

rule metaspades_assembly_from_unmapped:
    input:
        'Assembly/intermediary/{type}.r1.unmapped.fq',
        'Assembly/intermediary/{type}.r2.unmapped.fq',
        'Assembly/intermediary/{type}.se.unmapped.fq'
    output:
        directory('Assembly/intermediary/{type}.metaspades_unmapped.2'),
        'Assembly/intermediary/{type}.metaspades_unmapped.2/contigs.fa',
        'Assembly/intermediary/{type}.metaspades_unmapped.2.fa'
    params:
        runtime = "120:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ROOT_DIR + "/envs/IMP_assembly.yaml"
    log: "logs/assembly_metaspades_assembly_from_unmapped.{type}.log"
    message: "metaspades_assembly_from_unmapped: Performing {wildcards.type} assembly step 2 from unmapped reads using METASPADES"
    shell:
        METASPADES_ASSEMBLY_SHELL
