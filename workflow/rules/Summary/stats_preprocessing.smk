input_PREPRO = []
input_PREPRO.append(expand("Stats/{type}/{type}.read_counts.txt", type=TYPES))
if MG:
    input_PREPRO.append(expand("Stats/mg/mg.{read}_nonpareil.diversity.txt", read=["r1","r2"]))


rule PREPROCESSING:
    input:
        input_PREPRO
    output:
        touch('status/preprocessing.analysed')
