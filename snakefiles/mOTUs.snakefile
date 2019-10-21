rule allmotu:
	input:
		expand("mOTUs/{type}.mOTU.{measure}.tsv", type=TYPES, measure=['abundances','counts'])
	output:
		"mOTUs.done"
	shell:
		"""touch {output}"""

rule motus:
	input:
		"Preprocessing/{type}.r1.preprocessed.fq",
		"Preprocessing/{type}.r2.preprocessed.fq"
	output:
		"mOTUs/{type}",
		"mOTUs/{type}/RESULTS/mOTU.abundances.gz",
		"mOTUs/{type}/RESULTS/mOTU.counts.gz"
	threads: 12
	params:
		"{type}.r1.preprocessed.fq",
		"{type}.r2.preprocessed.fq"
	shell:
		"""mkdir -p {output[0]}
		cd {output[0]}
		ln -s ../../{input[0]}
		ln -s ../../{input[1]}
		mOTUs.pl --processors={threads} {params[0]} {params[1]}
		cd ../../"""

rule getResults:
	input:
		"mOTUs/{type}/RESULTS/mOTU.{measure}.gz"
	output:
		"mOTUs/{type}.mOTU.{measure}.tsv"
	shell:
		"""zcat {input} >> {output}"""

