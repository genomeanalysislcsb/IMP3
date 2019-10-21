rule proteomicsvarprep:
	input:
		"Analysis/annotation/proteinPredictions.mgmt.variants.human.formatted.faa"
	output:
		"proteomicsPreparation.done"
	shell:
		"""touch {output}"""

rule proteomicsvarformatting:
	input:
		"Analysis/annotation/proteinPredictions.mgmt.variants.human.faa"
	output:
		"Analysis/annotation/proteinPredictions.mgmt.variants.human.formatted.faa"
	shell:
		"""sed -s 's/[*]//' {input} | fasta_formatter -w 60 >> {output}"""


rule mergeProteomicsvarHuman:
	input:
		"Analysis/annotation/mgmt.assembly.merged.mgmt.variants.proteins.faa"
	output:
		"Analysis/annotation/proteinPredictions.mgmt.variants.human.faa"
	shell:
		"""cat {input} {config[humanProteome]} >> {output}"""
		
		
rule insertVariants:
	input:
		annotation="Analysis/annotation/annotation_CDS_RNA_hmms.gff",
		assembly="Assembly/mgmt.assembly.merged.fa",
		mgvar="Analysis/mg.variants.samtools.vcf.gz",
		mtvar="Analysis/mt.variants.samtools.vcf.gz",
	output:
		predictions="Analysis/annotation/mgmt.assembly.merged.mgmt.variants.filtered.proteins.faa",
		peptideTab="Analysis/annotation/mgmt.assembly.merged.mgmt.variants.filtered.pos.tsv",
		addInfo="Analysis/annotation/mgmt.assembly.merged.mgmt.variants.filtered.info"
	log:
		"Analysis/proteomics.log"
	shell:
		"""variant_monsterProtein_anvioFilter_filter.pl -a {input.annotation} -s {input.assembly} -v {input.mgvar} -v {input.mtvar} -n Analysis/annotation/mgmt.assembly.merged.mgmt.variants.filtered -p -anvio --min 7 >> {log} 2>&1"""

 
