METASPADES_ASSEMBLY_SHELL = """
if [ -d "{output[0]}" ]; then
    rm -rf {output[0]}
fi
rnaspades.py --meta \
 --pe1-1 {input[0]} \
 --pe1-2 {input[1]} \
 --pe1-s {input[2]} \
  -t {threads} \
  -m {BIGMEMTOTAL} \
  -k {KMER_STEPS} \
  {STRANDED_ARG} \
  -o {output[0]} > {log} 2>&1
ln -fs  {output[1]} {output[2]}
"""

kmersteps=range(config['assembly']['mink'],config['assembly']['maxk']+1,config['assembly']['step'])
KMER_STEPS =",".join(map(str,kmersteps))
STRANDED_ARG = ""
if MT_STRANDED == 1:
    STRANDED_ARG =  "-ss-rf"
elif MT_STRANDED == 2:
    STRANDED_ARG =  "-ss-fr" 

rule metaspades_assembly_from_preprocessing:
    input:
        'Preprocessing/mt.r1.preprocessed.fq',
        'Preprocessing/mt.r2.preprocessed.fq',
        'Preprocessing/mt.se.preprocessed.fq'
    output:
        directory('Assembly/intermediary/mt.metaspades_preprocessed.1'),
        'Assembly/intermediary/mt.metaspades_preprocessed.1/transcripts.fasta',
        'Assembly/intermediary/mt.metaspades_preprocessed.1.fa'
    resources:
        runtime = "48:00:00",
        mem_mb = BIGMEMCORE_MB
    threads: getThreads(BIGCORENO, bigmem=True)
    conda: ENVDIR + "/IMP_assembly.yaml"
    log: "logs/assembly_metaspades_assembly_from_preprocessing.mt.log"
    message: "metaspades_assembly_from_preprocessing: Performing mt assembly step 1 from preprocessed reads using MetaSpades"
    shell:
        METASPADES_ASSEMBLY_SHELL

rule metaspades_assembly_from_unmapped:
    input:
        'Assembly/intermediary/mt.r1.unmapped.fq',
        'Assembly/intermediary/mt.r2.unmapped.fq',
        'Assembly/intermediary/mt.se.unmapped.fq'
    output:
        directory('Assembly/intermediary/mt.metaspades_unmapped.2'),
        'Assembly/intermediary/mt.metaspades_unmapped.2/transcripts.fasta',
        'Assembly/intermediary/mt.metaspades_unmapped.2.fa'
    resources:
        runtime = "48:00:00",
        mem_mb = BIGMEMCORE_MB
    threads: getThreads(BIGCORENO, bigmem=True)
    conda: ENVDIR + "/IMP_assembly.yaml"
    log: "logs/assembly_metaspades_assembly_from_unmapped.mt.log"
    message: "metaspades_assembly_from_unmapped: Performing mt assembly step 2 from unmapped reads using METASPADES"
    shell:
        METASPADES_ASSEMBLY_SHELL
