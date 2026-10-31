
rule no_filtering:
    input:
        no_filtering_input_SE
    output:
        'Preprocessing/{type}.r1.preprocessed.fq.gz',
        'Preprocessing/{type}.r2.preprocessed.fq.gz',
        'Preprocessing/{type}.se.preprocessed.fq.gz'
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    # log:
    #     "logs/preprocessing/no_filtering.{type}.log"
    message:
        "no_filtering: not filtering {wildcards.type} reads."
    shadow: "copy-minimal"
    shell:
        """
        R1=mg.r1.fq
        R2=mg.r2.fq
        SE=`echo {input[0]} | sed 's/Preprocessing\///g'`

        L_R1=`echo {output[0]} | sed 's/Preprocessing\///g'`
        L_R2=`echo {output[1]} | sed 's/Preprocessing\///g'`
        L_SE=`echo {output[2]} | sed 's/Preprocessing\///g'`

        cd Preprocessing
        ln -fs ${{R1}} ${{L_R1}} && touch -h ${{L_R1}}
        ln -fs ${{R2}} ${{L_R2}} && touch -h ${{L_R2}}
        ln -fs ${{SE}} ${{L_SE}} && touch -h ${{L_SE}}
        cd ..
        """
