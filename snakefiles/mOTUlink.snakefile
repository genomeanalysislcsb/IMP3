rule allmotulinks:
	input:
		expand("mOTU_links/{COG}.bestPresentHitPhylogeny.tsv",COG=config["COGS"])
	output:
		"mOTU_links/bestPresentHitPhylogeny.tsv"
	shell:
		"""
		cat <(head -n 1 {input} | tail -n 1) <(cat {input} | grep -v queryID ) > {output}
		"""

rule fetchmg:
	input:
		"Analysis/annotation/genepredictions.Prodigal/prodigal.C.fna",
		"Analysis/annotation/prodigal.C.faa",
		"mOTU_links/"
	output:
		"mOTU_links/marker_genes",
		expand("mOTU_links/marker_genes/{COG}.fna",COG=config["COGS"])
	threads: 12
	shell:
		"""
		fetchMG.pl -o {output[0]} -t {threads} -d {input[0]}  -m extraction {input[1]}
		"""

	
rule get_present_mOTUs:
	input:
		expand("mOTUs/{type}.mOTU.counts.tsv",type=TYPES)
	output:
		"mOTU_links/",		
		"mOTU_links/all.mOTUs.present"
	shell:
		"""
	#	mkdir {output[0]}
		grep "motus.processing" -v {input} -h | cut -f 1 | sort | uniq > {output[1]}
		
		"""

rule extract_mOTU_db_genes:
	input:
		"mOTU_links/all.mOTUs.present"
	output:
		"mOTU_links/{COG}.present.mOTU.DBgenenames.txt",
		"mOTU_links/{COG}.present.DBgenes.fna"
	params:
		"{COG}"
	shell:
		"""
		Rscript $imp_root/src/mergeMOTUgeneInfo.R {params} {input} {output[0]}
		fastaExtract.pl $imp_lib/motus_linking/{params}/{params}.mOTU.v1.padded {output[0]} > {output[1]}
		"""

rule blast_mOTU_db_genes:
	input:
		"mOTU_links/{COG}.present.DBgenes.fna",
		"mOTU_links/marker_genes/{COG}.fna"
	output:
		"mOTU_links/{COG}.bestPresentHits.tsv",
		"mOTU_links/{COG}.bestPresentHitPhylogeny.tsv"
	params:
		"{COG}"
	shell:
		"""
		makeblastdb -in {input[0]} -dbtype nucl -parse_seqids
		blastn -db {input[0]} -query {input[1]} -outfmt 7 -out {output[0]} -max_target_seqs 1
		Rscript $imp_root/src/mergeMOTUblast.R {params} {output[0]} {output[1]}
		"""


