rule maxbin_get_depth:
    input:
        "Stats/mg/%s.assembly.contig_depth.txt" % ASS,
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Stats/mg/%s.assembly.contig_depth.0.txt" % ASS
    threads: 1
    params:
        runtime = "3:00:00",
        mem = MEMCORE
#    log: "logs/binning_maxbin_get_depth.log"
    message: "maxbin_get_depth: Getting depth in Maxbin format."
    shell:
        """
        cat {input[0]} \
         <(awk 'BEGIN {{OFS="\t"}}; {{print $0, 0}}' \
         <(diff --new-line-format="" --unchanged-line-format="" \
         <(grep "^>" {input[1]} |sed -e 's/>//g' | sort) \
         <(cut -f1 {input[0]} | sort))) > {output[0]}
        """
