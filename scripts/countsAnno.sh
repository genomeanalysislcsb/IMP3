#! /bin/bash -l

featureCount=/home/users/aheintzbuschart/bin/subread-1.5.1-Linux-x86_64/bin/featureCounts
threads=12
annAdd="/home/users/aheintzbuschart/myScripts/hmmscan_addBest2gff.pl"
gff=annotation.filt.gff
result=annotation.filt
fasta=prokka.faa
cd Analysis/annotation

perl $annAdd -h prokka.faa.kegg.hmmscan -a $gff -n KEGG -o $result.K.gff -g $(grep ">" $fasta | wc -l) -t
perl $annAdd -h prokka.faa.metacyc.hmmscan -a $result.K.gff -n metaCyc -o $result.K.M.gff -g $(grep ">" $fasta | wc -l) -t
perl $annAdd -h prokka.faa.Pfam-A.hmmscan -a $result.K.M.gff -n Pfam -o $result.K.M.P.gff -g $(grep ">" $fasta | wc -l)
perl $annAdd -h prokka.faa.tigrpfam.hmmscan -a $result.K.M.P.gff -n TIGR -o $result.K.M.P.T.gff -g $(grep ">" $fasta | wc -l) -t
perl $annAdd -h prokka.faa.swissprot.hmmscan -a $result.K.M.P.T.gff -n swissProt -o $result.allHMM.gff -g $(grep ">" $fasta | wc -l) -t

for ome in mg mt
do
$featureCount -p -O -t CDS -g KEGG -o $ome.KEGG.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
$featureCount -p -O -t CDS -g metaCyc -o $ome.metaCyc.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
$featureCount -p -O -t CDS -g Pfam -o $ome.Pfam.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
$featureCount -p -O -t CDS -g TIGR -o $ome.TIGR.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
$featureCount -p -O -t CDS -g swissProt -o $ome.swissProt.counts.txt -a annotation.filt.allHMM.gff -T $threads ../../Assembly/$ome.reads.sorted.bam
done

