rule proteomicsprep:
	input:
		"Analysis/annotation/prodigal.noEnds.human.faa"
	output:
		"proteomics.prepped"
	shell:
		"""touch {output}"""

rule cutends:
	input:
		"Analysis/annotation/prodigal.C.faa"
	output:
		"Analysis/annotation/prodigal.noEnds.faa"
	shell:
		"""trypsinStartEndFastaOnly.filter2.pl {input} | fasta_formatter -w 60 >> {output}"""


rule mergeHuman:
	input:
		"Analysis/annotation/prodigal.noEnds.faa"
	output:
		"Analysis/annotation/prodigal.noEnds.human.faa"
	shell:
		"""cat {input} {config[imp_lib]}/human.cds.corrected.protein.fa >> {output}"""
 
