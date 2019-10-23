rule binny:
	input:
		expand("Analysis/results/{ASS}_results.Rdat", ASS=ASS),
		expand("Binning/{ASS}.vizbin.with-contig-names.points", ASS=ASS),
		"Analysis/annotation/annotation.CDS.RNA.essential.gff"
	output:
		expand("Binning/reachabilityDistanceEstimates_contigs.{pk}.{nn}.tsv \
		Binning/contigs2clusters.{pk}.{nn}.tsv \
		Binning/contigs2clusters.{pk}.{nn}.RDS \
		Binning/clusteringWS_contigs.{pk}.{nn}.Rdata \
		Binning/finalClusterMap_contigs.{pk}.{nn}.pdf".split(),pk=config["binning"]["pk"],nn=config["binning"]["nn"])
	shell:
		"""
		mkdir -p Binning/cluster.{config[binning][pk]}.{config[binning][nn]}.Files
		Rscript $imp_root/src/binny.R {input} Binning {config[binning][pk]} {config[binning][nn]}
		"""

rule cut_rRNA:
	input:
		"Assembly/{ASS}.assembly.merged.fa",
		"Analysis/annotation/annotation.filt.gff"
	output:
		"Binning/{ASS}.assembly.merged.cut.fa"
	params:
		"Binning/rRNAcutting.log"
	shell:
		"""fastaExtractCutRibosomal1000.pl -f {input[0]} -g {input[1]} -l {params} -o {output} -c {config[vizbin][cutoff]}"""

rule rerun_vizbin:
	input:
		"Binning/{ASS}.assembly.merged.cut.fa"
	output:
		"Binning/{ASS}.vizbin.with-contig-names.points"
	shell:
		"""
		echo "[x] VIZBIN `date +"%Y/%m/%d %H:%M:%S"`"
		TMP_VIZBIN=$(mktemp -dt "VIZBIN_XXXXXX")
		java -jar {config[vizbin][jarfile]} \
		-a {config[vizbin][dimension]} \
		-c {config[vizbin][cutoff]} \
		-i {input} \
		-o $TMP_VIZBIN/data.points \
		-k {config[vizbin][kmer]} \
		-p {config[vizbin][perp]}

		if [ -f $TMP_VIZBIN/data.points ]
            then
                paste <(grep "^>" {input} | sed -e 's/>//') \
                <(cat $TMP_VIZBIN/data.points | sed -e 's/,/\t/') > {output}
            fi
        rm -rf $TMP_VIZBIN
		"""
