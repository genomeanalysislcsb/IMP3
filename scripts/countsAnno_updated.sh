#! /bin/bash -l
#this script should work with the output of the output of the annotate.snakefile
#if you've run Susana's version of that workflow, use annotation_CDS_RNA_hmms.gff for $annotation_file 

featureCount=/work/projects/ecosystem_biology/local_tools/subread-1.5.1-Linux-x86_64/bin/featureCounts
threads=12
annotation_file=annotation.CDS.RNA.hmms.gff
#if you've run Susana's version of that workflow, use annotation_CDS_RNA_hmms.gff for $annotation_file 

#annAdd="/work/projects/ecosystem_biology/local_tools/IMPPatrick/IMPnew/bin/hmmscan_addBest2gff.pl"
#gff=annotation.filt.gff
#result=annotation.filt
#fasta=prokka.faa
cd Analysis/annotation

#perl $annAdd -h prokka.faa.kegg.hmmscan -a $gff -n KEGG -o $result.K.gff -g $(grep ">" $fasta | wc -l) -t
#perl $annAdd -h prokka.faa.metacyc.hmmscan -a $result.K.gff -n metaCyc -o $result.K.M.gff -g $(grep ">" $fasta | wc -l) -t
#perl $annAdd -h prokka.faa.Pfam-A.hmmscan -a $result.K.M.gff -n Pfam -o $result.K.M.P.gff -g $(grep ">" $fasta | wc -l)
#perl $annAdd -h prokka.faa.tigrpfam.hmmscan -a $result.K.M.P.gff -n TIGR -o $result.K.M.P.T.gff -g $(grep ">" $fasta | wc -l) -t
#perl $annAdd -h prokka.faa.swissprot.hmmscan -a $result.K.M.P.T.gff -n swissProt -o $result.allHMM.gff -g $(grep ">" $fasta | wc -l) -t

for ome in mg mt
do
for db in KEGG metaCyc Pfam-A TIGR swissProt Cas essential
do
db_id=${db//-/_}
$featureCount -p -O -t CDS -g $db_id -o $ome.$db.counts.txt -a $annotation_file -T $threads ../../Assembly/$ome.reads.sorted.bam
#you could also use -s 2 for stranded mt

#$featureCount -p -O -t CDS -g metaCyc -o $ome.metaCyc.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
#$featureCount -p -O -t CDS -g Pfam -o $ome.Pfam.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
#$featureCount -p -O -t CDS -g TIGR -o $ome.TIGR.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
#$featureCount -p -O -t CDS -g swissProt -o $ome.swissProt.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
done
done

