#!/bin/bash -l

input0="Binning/selected_DASTool_bins"
output="Binning/selected_genome_locations.tsv"
      
for cf in `ls $input0 | grep -v contigs.fa`
do
	printf "%s\t%s\n" "Binning/selected_DASTool_bins/$cf/contigs.fa" "$cf"
done >$output

input="Binning/selected_genome_locations.tsv"
output0="Binning/selected_GTDBtk.done"
output1="Binning/selected_GTDBtk_classifications"
threads=8

__conda_setup="$('/work/projects/ecosystem_biology/local_tools/miniconda/miniconda2/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"

conda activate gtdbtk
mkdir -p  $output1

for genome in `cut -f 2 $input`
do
	gdirp=`grep $genome -w $input | cut -f 1`
        gdir=${gdirp/contigs.fa/}
	echo "# $genome gtdbtk classify start"       
	gtdbtk classify_wf --genome_dir $gdir -x fa --out_dir ${gdir}GTDB --cpus ${threads}

	if [ -f ${gdir}GTDB/gtdbtk.bac120.summary.tsv ]
        then
        	tail -n 1 ${gdir}GTDB/gtdbtk.bac120.summary.tsv | sed "s/contigs/$genome/" >> ${output1}/gtdbtk.bac120.summary.tsv
        fi
        if [ -f ${gdir}GTDB/gtdbtk.ar122.summary.tsv ]
        then
        	tail -n 1 ${gdir}GTDB/gtdbtk.ar122.summary.tsv | sed "s/contigs/$genome/" >> ${output1}/gtdbtk.ar122.summary.tsv
        fi
	echo "# $genome gtdbtk classify done" 
done
touch ${output0}
conda deactivate

