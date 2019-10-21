rule rerunProdigal:
	input:
		expand("Assembly/{ass}.assembly.merged.fa", ass=ASS)
	output:
		genes="Analysis/annotation/prodigal.C.faa",
		gff="Analysis/annotation/prodigal.C.gff",
		prodigalout="Analysis/annotation/genepredictions.Prodigal"
	log:
		"Analysis/annotation/genepredictions.Prodigal/gene.prediction.log"
	shell:
		"""prodigal -f gff -a {output.prodigalout}/prodigal.C.faa.tmp -d {output.prodigalout}/prodigal.C.fna.tmp -p meta -o {output.prodigalout}/prodigal.C.gff -i {input} -q >> {log} 2>&1
		MOCATGenePredictionProdigal_fixed_aux.pl {output.prodigalout}/prodigal.C.fna.tmp {output.prodigalout}/prodigal.C.faa.tmp {output.prodigalout}/prodigal.C.tab >> {log} 2>&1 
		mv {output.prodigalout}/prodigal.C.faa {output.genes}
		mv {output.prodigalout}/prodigal.C.gff {output.gff}"""

