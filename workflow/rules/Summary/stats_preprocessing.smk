rule PREPROCESSING:
    input:
        expand("Stats/{type}/{type}.read_counts.txt", type=TYPES)
    output:
        touch('status/preprocessing.analysed')
