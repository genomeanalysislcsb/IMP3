rule call_gene_depth:
    input:
        "Analysis/annotation/annotation.filt.gff",
        "Assembly/{type}.reads.sorted.bam"
    output:
        "Stats/{type}/annotation/{type}.annotation_cov.bed",
        "Stats/{type}/annotation/{type}.gene_depth.hist",
        "Stats/{type}/annotation/{type}.gene_depth.avg",
        "Stats/{type}/annotation/{type}.gene.len"
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/analysis_call_gene_depth.{type}.log"
    message: "call_gene_depth: Getting data on gene coverage with {wildcards.type} reads."
    shell:
        """
        coverageBed -hist -abam {input[1]} -b <(awk '$4 < $5' {input[0]}) | grep -v "^all" > {output[0]}
        paste <(cat {output[0]} | cut -f9 | cut -f1 -d \";\" | sed -e \"s/ID=//g\") \
         <(cut -f10,11,12,13 {output[0]}) > {output[1]}
        ## This code was adapted and modified from the CONCOCT script to calculate depth
        ## It prints out a file that contains the average depth of all the genes
        awk -v OFS='\t' 'BEGIN {{pc=""}}
        {{
            c=$1;
            if (c == pc) {{
                    cov=cov+$2*$5;
            }} else {{
                print pc,cov;
                cov=$2*$5;
                pc=c
            }}
        }} END {{print pc,cov}}' < {output[1]} | tail -n +2 > {output[2]}


        # Record gene length file
        cut -f 1,4 {output[1]} | uniq > {output[3]}
        """
