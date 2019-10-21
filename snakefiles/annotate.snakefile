
rule mergegff:
	input:
		expand("Analysis/annotation/annotation.CDS.RNA.{db}.gff",db=config["hmmdbs"])
	output:
		"Analysis/annotation/annotation.CDS.RNA.hmms.gff"
	shell:
		"""mergegffs.pl {output} {input}"""

rule makegff:
	input:
		"Analysis/annotation/prodigal.C.faa.{db}.hmmscan",
		"Analysis/annotation/annotation.CDS.RNA.gff",
		"Analysis/annotation/prodigal.C.faa"
	output:
		"Analysis/annotation/annotation.CDS.RNA.{db}.gff"
	params:
		"{db}",
		lambda wildcards: config["hmmdbs"][wildcards.db]["trim"]
	shell:
		"""hmmscan_addBest2gff.pl -i {input[0]} -a {input[1]} -n {params[0]} {params[1]} -o {output} -g $(grep ">" {input[2]} | wc -l)"""

rule hmmer:
	input:
		"Analysis/annotation/prodigal.C.faa"
	output:
		"Analysis/annotation/prodigal.C.faa.{db}.hmmscan"
	threads: 12
	params: 
		"/work/projects/ecosystem_biology/data/hmm/{db}",
		lambda wildcards: config["hmmdbs"][wildcards.db]["cutoff"]
	shell:
		"""hmmsearch --cpu {threads} {params[1]} --noali --notextw --tblout {output} {params[0]}/*.hmm {input} >/dev/null"""

rule mergeGenesRNA:
    input:
        "Analysis/annotation/annotation.filt.gff",
        "Analysis/annotation/prodigal.C.gff"
    output:
        "Analysis/annotation/annotation.CDS.RNA.gff"
    shell:
        """{{ grep Prodigal -v {input[0]}; grep "^#" -v {input[1]} | sed 's/ID=/ID=contig_gene/'; }} > {output}"""

