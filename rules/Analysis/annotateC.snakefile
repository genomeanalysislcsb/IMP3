rule hmmer:
	input:
		"Analysis/annotation/prokka.faa"
	output:
		"Analysis/annotation/prokka.faa.{db}.hmmscan"
#	threads: 12
	params: 
		"/work/projects/ecosystem_biology/data/hmm/{db}",
		lambda wildcards: config["hmmdbs"][wildcards.db]["cutoff"]
	shell:
		"""hmmsearch --cpu {threads} {params[1]} --noali --notextw --tblout {output} {params[0]}/*.hmm {input} >/dev/null"""

rule makegff:
        input:
                "Analysis/annotation/prokka.faa.{db}.hmmscan",
                "Analysis/annotation/annotation.filt.gff",
                "Analysis/annotation/prokka.faa"
        output:
                "Analysis/annotation/annotation.CDS.RNA.{db}.gff"
        params:
                "{db}",
                lambda wildcards: config["hmmdbs"][wildcards.db]["trim"]
        shell:
                """hmmscan_addBest2gff.pl -i {input[0]} -a {input[1]} -n {params[0]} {params[1]} -o {output} -g $(grep ">" {input[2]} | wc -l)"""

rule mergegff:
        input:
                expand("Analysis/annotation/annotation.CDS.RNA.{db}.gff",db=config["hmmdbs"])
        output:
                "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
        shell:
                """mergegffs.pl {output} {input}"""

