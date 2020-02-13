rule megahit_hybrid_assembly_1:
    input:
        'Preprocessing/mg.r1.preprocessed.fq',
        'Preprocessing/mg.r2.preprocessed.fq',
        'Preprocessing/mg.se.preprocessed.fq',
        'Preprocessing/mt.r1.preprocessed.fq',
        'Preprocessing/mt.r2.preprocessed.fq',
        'Preprocessing/mt.se.preprocessed.fq',
        'Assembly/intermediary/mt.megahit_preprocessed.1/final.contigs.fa',
        'Assembly/intermediary/mt.megahit_unmapped.2/final.contigs.fa',
    output:
        'Assembly/intermediary/mgmt.megahit_hybrid.1/final.contigs.fa',
        'Assembly/intermediary/mgmt.megahit_hybrid.1.fa',
        directory('Assembly/intermediary/mgmt.megahit_hybrid.1')
    params:
        runtime = "120:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ENVDIR + "/IMP_assembly.yaml"
    log: "logs/assembly_megahit_hybrid_assembly_1.log"
    message: "megahit_hybrid_assembly_1: Performing hybrid assembly 1 from preprocessed reads using MEGAHIT"
    shell:
        """
        if [ -d "{output[2]}" ]; then
            rm -rf {output[2]}
        fi
        MAX_MEM="$(({BIGMEMTOTAL} * 1000000000))"
        megahit -1 {input[0]},{input[3]} \
         -2 {input[1]},{input[4]} \
         -r {input[2]},{input[5]},{input[6]},{input[7]} \
         -o {output[2]} \
         --k-min {config[assembly][mink]} \
         --k-max {config[assembly][maxk]} \
         --k-step {config[assembly][step]} \
         -t {threads} --cpu-only \
         -m "${{MAX_MEM}}" \
         --mem-flag 1 > {log} 2>&1
        ln -fs $(echo {output[0]} | cut -f 3,4 -d /) {output[1]} && touch -h {output[1]}
        """

rule megahit_hybrid_assembly_2:
    input:
        'Assembly/intermediary/mgmt.r1.megahit_hybrid.mt.unmapped.fq',
        'Assembly/intermediary/mgmt.r2.megahit_hybrid.mt.unmapped.fq',
        'Assembly/intermediary/mgmt.se.megahit_hybrid.mt.unmapped.fq',
        'Assembly/intermediary/mgmt.r1.megahit_hybrid.mg.unmapped.fq',
        'Assembly/intermediary/mgmt.r2.megahit_hybrid.mg.unmapped.fq',
        'Assembly/intermediary/mgmt.se.megahit_hybrid.mg.unmapped.fq',
    output:
        'Assembly/intermediary/mgmt.megahit_hybrid.2/final.contigs.fa',
        'Assembly/intermediary/mgmt.megahit_hybrid.2.fa',
        directory('Assembly/intermediary/mgmt.megahit_hybrid.2')
    params:
        runtime = "120:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ENVDIR + "/IMP_assembly.yaml"
    log: "logs/assembly_megahit_hybrid_assembly_2.log"
    message: "megahit_hybrid_assembly_2: Performing hybrid assembly 2 from preprocessed reads using MEGAHIT"
    shell:
        """
        if [ -d "{output[2]}" ]; then
            rm -rf {output[2]}
        fi

        MAX_MEM="$(({BIGMEMTOTAL} * 1000000000))"

        echo "Performing second hyrbid assembly step using MEGAHIT"
        megahit -1 {input[0]},{input[3]} \
         -2 {input[1]},{input[4]} \
         -r {input[2]},{input[5]} \
         -o {output[2]} \
         --k-min {config[assembly][mink]} \
         --k-max {config[assembly][maxk]} \
         --k-step {config[assembly][step]} \
         -t {threads} --cpu-only \
         -m ${{MAX_MEM}} \
         --mem-flag 1 > {log} 2>&1
        ln -fs $(echo {output[0]} | cut -f 3,4 -d /) {output[1]} && touch -h {output[1]}
        """
