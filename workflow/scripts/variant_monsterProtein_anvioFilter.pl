#!/usr/bin/env perl -w

# script to include variant proteins in protein database
# takes 6 inputs: -a, a table with the gene positions, starts and stops (like the filtered .gff output of IMP);
#		  -s, a fasta file with the contigs
#                 -v, a list of .vcf files of the same contigs
#		  -p, if not given, no fasta files will be generated
#		  -min, minimum peptide length for the additional peptides and to decide which proteins to keep, defaults to 5
#		  -max, maximum peptide length for the additional peptides, defaults to 144
#                 -anvio, if set, vcf files with the word samtools in the file name will be filtered according to anvio heuristics (at least the first), otherwise the phred score is used to 99% certainty
#                 -varmax, the number of variants to consider per peptide, default is 60000
#			-u the minimum number of uniqpeptides
# 3 outputs: - a fasta file with the amino acid sequences of all peptides, including variants; <name>.peptides.variants.faa
#	     - a table with the peptides and variant positions; for original peptides, one position without change is indicated
#	     - a log file <name>.variants.log
#

#written by Anna Heintz-Buschart (Feb 2016) based on code by Nicolas Pinel (July 04, 2013)


use strict;
use Getopt::Long;
use Bio::SeqIO;

# Obtain inputs
my ($ann,$seq,@vcf,$print,$uniqpeps,$minPeplength,$maxPeplength,$maxvar,$anvioheur,$nameout);
GetOptions('a|ann=s' => \$ann,
	   's|seq=s' => \$seq,
	   'v|vcf=s' => \@vcf,
	   'o|outname=s' => \$nameout,
	   'p|print' => \$print,
	   'u|uniqpeps=i' => \$uniqpeps,
	   'min=i' => \$minPeplength,
	   'max=i' => \$maxPeplength,
	   'varmax=i' => \$maxvar,
	   'anvio' => \$anvioheur,);
$minPeplength = $minPeplength ? $minPeplength : 5;
$maxPeplength = $maxPeplength ? $maxPeplength : 144;
$nameout = $nameout ? $nameout : "variantpredictions";
$uniqpeps = $uniqpeps ? $uniqpeps : 2;
$maxvar = $maxvar ? $maxvar : 60000;

# Load sequences
my $seqbank = &loadSequences($seq);

# Load annotation (tab file)
my $annot = &loadAnnotation($ann);

# Load variant file (subroutine written)
foreach my $vcf (@vcf){
    if (index($vcf, "samtools") != -1) {
	if ($anvioheur) {
	    &loadVariantsAnvST($vcf);
	} else{
	    &loadVariantsST($vcf);
	}
    } else{
	&loadVariants($vcf);
    }
}

# Load amino acid table (manually written subroutine)
my $aa_table = &loadAAreg();
my $aab_table = &loadAA();
my $fix = 0;
my $prov_out = $nameout.".proteins.faa";
#(my $pepv_out = $seq) =~ s/\.\w+$/\.$nameout.addpeptides.faa/;
my $log_out = $nameout.".info";
my $tab_out = $nameout.".pos.tsv";

# Print output files
my ($o_prov, $o_pepv, $o_log, $o_tab);
if ($print) {
    open $o_prov, '>', $prov_out;
    #open $o_pepv, '>', $pepv_out;
    open $o_log, '>', $log_out;
    open $o_tab, '>', $tab_out;
}
foreach my $k (sort keys %{$annot}) { #for every contig
    foreach my $feat (@{$$annot{$k}}) { #for every gene
		my $id = $$feat{'gene'}; #gene name
		if ($$feat{'variants'}) { #variants in gene
			my ($prov,$idv);
			my $nuc = &getSeqProtein($k, $feat); #sequence retrieval with contig and gene info, gets CDS
			my @pro = &translateProtein($nuc, $k, $feat); #get translated protein, includes removing ends for incomplete
			my $proseq = shift @pro; #the first element is the protein
			print $o_prov qq(\>$id\n$proseq\n) if ($print && $proseq ne ""); #print translated protein to output fasta
			my @varcoll = &getAllVars($feat,$$feat{'variants'}); #get variants within protein - sorted
			my $vnum = $#varcoll; 
			my $vcnt = 0;
			my $pcnt = 0;
			my $shiftPro = 0;
			my $pep1 = 1; #set to 1 when it's the first peptide to be printed for a protein
			my %peptidesForOut = ();
			my @fscoll = ();
			## and now print variants by possibilites per feature
			while ($vcnt <= $vnum){ #for all variant sets
				my $pv = 0;
				my %peptides = &getPeptides($k, $feat, $varcoll[$vcnt]); #get all peptides, based on the variants
				foreach my $pep (keys %peptides) { #for every peptide from the variants
					my $pepex = 0;
					foreach my $opep (@pro) { #check if peptide is a peptide of the original protein (the rest of pro contains the peptides)
						if($opep eq $pep){
							$pepex = 1; #if it is, nothing is done
							last;
						}
					}
					if ($pepex == 0 && length($pep)>=$minPeplength && length($pep)<=$maxPeplength){ #here we check for length
						if (!$peptidesForOut{$pep}) { #check that the peptide isn't already in the variants of this protein
							$peptidesForOut{$pep} = 1; #list of variant peptides in this protein 
							my $pvi = 0;
							my %outinf = ();
							delete $peptidesForOut{$pep} unless ($peptides{$pep}[$pvi]); #only keep peptide if there's all the data for this peptide
							while ($peptides{$pep}[$pvi]){
								$outinf{$peptides{$pep}[$pvi]{'origpos'}."\t".$peptides{$pep}[$pvi]{'origref'}."\t".$peptides{$pep}[$pvi]{'orignew'}} = 1;
								$pvi++;
							}
							foreach my $outf (keys %outinf){
								print $o_tab qq($pep\t$id\t$k\t$outf\tvariant\n) if ($print); # the log
							}
						}
					}
				}
				$vcnt++;
			}
			if (keys(%peptidesForOut)) { #if we really have peptides
				$pcnt++;
				my $varcnt = 1;
				my @stoppeps = ();
				my $stpcnt = 0;
				my $normalpeps = "";
				### still to do: count peptides in the proteins, if we use two peptides for identification, we don't need to write single peptides
				if ($proseq eq "") { #in case there was no protein, but we have a peptide now
					my $npepcnt = 0;
					foreach my $pepv (keys %peptidesForOut) {
						if($pepv =~ /B/){
							$pepv =~ s/B/\*/g;
							push @stoppeps,$pepv;
							$stpcnt++;
						}else{
							$normalpeps .= $pepv;
							$npepcnt++;
						}
					}
					my $stps = $stpcnt >= 1 ? 1 : 0;
					if($npepcnt + $stps >= $uniqpeps){#only print if we have enough peptides
						while($stpcnt >0){
							print $o_prov qq(\>$id.varpep.$varcnt\n$normalpeps$stoppeps[$stpcnt-1]\n) if ($print);
							$stpcnt--;
							$varcnt++;
						}
					}
				}elsif(substr($proseq,-1) eq "*"){
					my $npepcnt = 0;
					foreach my $pepo (@pro) {
						if(substr($pepo,-1) eq "B" && length($pepo) >= $minPeplength){
							$pepo =~ s/B/\*/g;
							push @stoppeps,$pepo;
							$stpcnt++;
						}elsif(length($pepo) >= $minPeplength){
							$normalpeps .= $pepo;
							$npepcnt++;
						}
					}
					foreach my $pepv (keys %peptidesForOut) {
						if($pepv =~ /B/){
							$pepv =~ s/B/\*/g;
							push @stoppeps,$pepv;
							$stpcnt++;
						}else{
							$normalpeps .= $pepv;
							$npepcnt++;
						}
					}
					my $stps = $stpcnt >= 1 ? 1 : 0;
					if($npepcnt + $stps >= $uniqpeps){#only print if we have enough peptides
						while($stpcnt >0){
							print $o_prov qq(\>$id.varpep.$varcnt\n$normalpeps$stoppeps[$stpcnt-1]\n) if ($print);
							$stpcnt--;
							$varcnt++;
						}
					}
				}else{
					my $npepcnt = $#pro + 1;
					$normalpeps = $proseq;
					foreach my $pepv (keys %peptidesForOut) {
						if($pepv =~ /B/){
							$pepv =~ s/B/\*/g;
							push @stoppeps,$pepv;
							$stpcnt++;
						}else{
							$normalpeps .= $pepv;
							$npepcnt++;
						}
					}
					my $stps = $stpcnt >= 1 ? 1 : 0;
					if($npepcnt + $stps >= $uniqpeps){#only print if we have enough peptides
						while($stpcnt >0){
							print $o_prov qq(\>$id.varpep.$varcnt\n$normalpeps$stoppeps[$stpcnt-1]\n) if ($print);
							$stpcnt--;
							$varcnt++;
						}
					}
				}
				
						#### they should go in the same file, but we put the header in front of the protein without end, if it has a stopcodon, else in front of that protein
						#### also, we need to check for a stop codon in the peptides.
						##### we keep the peptides with stop-codon separate,
						#### then, we make a string with all the other peptides (if there are any), we add those to the original peptides (except with stop),
						#### add the peptides with the stop to the end (if there were any other additional peptides, also add the orginal one)
						### we don't need to check for P at the beginning of the peptides, because that can't happen
			}
		} else { #if there are no variants in the protein
			++$fix; #for the log
			my $nuc = &getSeqProtein($k, $feat); #just get the sequence
			my @pro = &translateProtein($nuc, $k,$feat); # translate
			my $proseq = shift @pro; #the first element is the protein
			print $o_prov qq(\>$id\n$proseq\n) if ($print && $proseq ne "" && $#pro >= $uniqpeps-1); #write, unless the sequence is empty or not enough peptides for identification
		}
    }
}
print qq(invariant\: $fix\n);
if ($print) {
    close $o_prov;
    #open $o_pepv, '>', $pepv_out;
    close $o_log;
    close $o_tab;
}
exit;

#### subroutines ###

sub usage {
    
    exit;
}

sub array_minus(\@\@) {
    my %e = map{ $_ => undef } @{$_[1]};
    return grep( ! exists( $e{$_} ), @{$_[0]} );
}

sub loadSequences {
    my $f = shift;
    my $inseq = Bio::SeqIO->new(-file => $f,
				-format => 'fasta');
    my %seqs = ();
    while (my $seqobj = $inseq->next_seq()) {
	$seqs{$seqobj->display_id()} = $seqobj->seq();
    }
    return \%seqs;
}

sub loadAnnotation {
    my $f = shift;
    open my $in, '<', $f;
    my @headers = ('gene','strand','start','stop','startcodon','stopcodon','completeness','startcodonvar','stopcodonvar','completenessvar');
    my %feats = ();
    while (<$in>) {
		last if (/##FASTA/);
		next if (/^\#/);
		chomp;
		my %hash = ();
		my @f = split(/\t/);
		my $ctg = $f[0];
		next if ($f[2] ne "CDS");
		my $attributes = $f[8];
		my @gene = split(';',$attributes);
		my $gene = "";
		my $partinf = "";
		my $startc = "";
		my $stopc = "";
		my $complinf = "";
		foreach my $afield (@gene){
			if($afield =~ /^ID=/){
				($gene = $afield) =~ s/^ID=//;
			}elsif($afield =~ /^partial=/) {
				($partinf = $afield) =~ s/^partial=//;
				if ($f[6] eq '-') {
					$partinf = join('', reverse(split(//, $partinf)));
				}
				if($partinf == "00"){
					$startc = "yes";
					$stopc = "yes";
					$complinf = "complete";
				}elsif($partinf == "01"){
					$startc = "yes";
					$stopc = "no";
					$complinf = "incomplete";
				}elsif($partinf == "10"){
					$startc = "no";
					$stopc = "yes";
					$complinf = "incomplete";
				}else{
					$startc = "no";
					$stopc = "no";
					$complinf = "incomplete";
				}
			}
		}
		#my $gene = $gene[0];
		#$gene =~ s/ID=//;
		#$ctg =~ s/_gene\d+$//g; #metagenemark naming
		
		@hash{@headers} = ($gene,$f[6],$f[3],$f[4],$startc,$stopc,$complinf,$startc,$stopc,$complinf);
		push(@{$feats{$ctg}}, \%hash);
    }
    return \%feats;
}

sub loadVariants {
    my $f = shift;
    if ($f =~ /gz$/) {
	(my $nf = $f) =~ s/\.gz$//;
	system("zcat $f > $nf");
	$f = $nf;
    }
    open my $in, '<', $f;
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my @f = split(/\t/);
	foreach my $feat (@{$$annot{$f[0]}}) {
	    if ($f[1] >= $$feat{'start'} && $f[1] <= $$feat{'stop'}) {
			my $c = $f[4] =~ tr/,//;
			if ($c >0) {
				my @allAlts = split(",",$f[4]);
				foreach my $alt (@allAlts){
					push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $alt});
				}
			} else{
				push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $f[4]});
			}
	    }
	}
    }
}

sub loadVariantsAnvST {
    my $f = shift;
    if ($f =~ /gz$/) {
	(my $nf = $f) =~ s/\.gz$//;
	system("zcat $f > $nf");
	$f = $nf;
    }
    open my $in, '<', $f;
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my @f = split(/\t/);
	my @info = split(";",$f[7]);
	my @DP4 = split("=",$info[2]);
	my @covs = split(",",$DP4[1]);
	my $nx = $covs[0] + $covs[1];
	my $ny = $covs[2] + $covs[3];
	my $totcov = $nx + $ny;
	my $n1 = $nx >= $ny ? $nx : $ny;
	my $n2 = $nx >= $ny ? $ny : $nx;
	my $nrat = $n2/$n1;
	my $heurexp = $totcov**(1/3);
	my $heurexp2 = $heurexp - 1.45;
	my $heur = (1/3)**$heurexp2;
	my $heurFin = $heur + 0.05;
	if ($nrat > $heurFin) {
	    foreach my $feat (@{$$annot{$f[0]}}) {
		if ($f[1] >= $$feat{'start'} && $f[1] <= $$feat{'stop'}) {
		    my $c = $f[4] =~ tr/,//;
		    if ($c >0) {
			my @allAlts = split(",",$f[4]);
			foreach my $alt (@allAlts){
			    push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $alt});
			}
		    } else{
			push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $f[4]});
		    }
		}
	    }
	}
    }
}

sub loadVariantsST {
    my $f = shift;
    if ($f =~ /gz$/) {
	(my $nf = $f) =~ s/\.gz$//;
	system("zcat $f > $nf");
	$f = $nf;
    }
    open my $in, '<', $f;
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my @f = split(/\t/);
	if ($f[5] >= 20) {
	    foreach my $feat (@{$$annot{$f[0]}}) {
		if ($f[1] >= $$feat{'start'} && $f[1] <= $$feat{'stop'}) {
		    my $c = $f[4] =~ tr/,//;
		    if ($c >0) {
			my @allAlts = split(",",$f[4]);
			foreach my $alt (@allAlts){
			    push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $alt});
			}
		    } else{
			push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $f[4]});
		    }
		}
	    }
	}
    }
}


sub loadAA {
    my %aa = ('ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'TTA' => 'L',
	      'TTG' => 'L', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TTT' => 'F', 'TTC' => 'F', 'ATG' => 'M',
	      'TGT' => 'C', 'TGC' => 'C', 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'GGT' => 'G', 'GGC' => 'G',
	      'GGA' => 'G', 'GGG' => 'G', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P', 'ACT' => 'T', 'ACC' => 'T',
	      'ACA' => 'T', 'ACG' => 'T', 'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	      'TAT' => 'Y', 'TAC' => 'Y', 'TGG' => 'W', 'CAA' => 'Q', 'CAG' => 'Q', 'AAT' => 'N', 'AAC' => 'N', 'CAT' => 'H',
	      'CAC' => 'H', 'GAA' => 'E', 'GAG' => 'E', 'GAT' => 'D', 'GAC' => 'D', 'AAA' => 'K', 'AAG' => 'K', 'CGT' => 'R',
	      'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TAA' => 'B', 'TAG' => 'B', 'TGA' => 'B');
    return \%aa;
}

sub loadAAreg {
    my %aareg = ('ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'TTA' => 'L',
	      'TTG' => 'L', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TTT' => 'F', 'TTC' => 'F', 'ATG' => 'M',
	      'TGT' => 'C', 'TGC' => 'C', 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'GGT' => 'G', 'GGC' => 'G',
	      'GGA' => 'G', 'GGG' => 'G', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P', 'ACT' => 'T', 'ACC' => 'T',
	      'ACA' => 'T', 'ACG' => 'T', 'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	      'TAT' => 'Y', 'TAC' => 'Y', 'TGG' => 'W', 'CAA' => 'Q', 'CAG' => 'Q', 'AAT' => 'N', 'AAC' => 'N', 'CAT' => 'H',
	      'CAC' => 'H', 'GAA' => 'E', 'GAG' => 'E', 'GAT' => 'D', 'GAC' => 'D', 'AAA' => 'K', 'AAG' => 'K', 'CGT' => 'R',
	      'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*');
    return \%aareg;
}



sub getSeqProtein {
    my ($ctg,$ref) = @_;
    my $length = $$ref{'stop'} - $$ref{'start'} + 1;
    my $seq = substr($$seqbank{$ctg}, $$ref{'start'} - 1, $length);
    $seq =~ tr/a-z/A-Z/;
    my $seqo = "";
    if ($$ref{'strand'} eq '-') {
		$seq = join('', reverse(split(//, $seq)));
		$seq =~ tr/ACGT/TGCA/;
		$seqo = $seq;
    } else {
		$seqo = $seq;
    }
    return $seqo;
}
 
 
sub translateProtein {
    my ($seq, $ctg, $ref) = @_;
    my $gene = $$ref{'gene'};
    $seq =~ tr/a-z/A-Z/;
    my @codons = $seq =~ /(.{3})/g;
    my $aa = '';
    foreach my $cod (@codons) {
		my $add = $$aa_table{$cod} ? $$aa_table{$cod} : 'X';
		$aa .= $add;
		if ($add eq "*") {
			last
		}
    }
    my $offset = 0;
#cut peptides at K or R unless P
	if ($$ref{'startcodon'} ne "yes" && $$ref{'stopcodon'} eq "yes") {
        my $ncut = 0;
        while ($aa ne "" && $ncut == 0) {
            my $k = index($aa, 'K') >= 0? index($aa, 'K') : length($aa)-1;
            my $r = index($aa, 'R') >= 0? index($aa, 'R') : length($aa)-1;
            my $cut = $k < $r ? $k : $r;
            $aa = substr($aa, $cut+1);
			$offset += 3*$cut+3;
			if ($aa eq "" || substr($aa,0,1) ne "P") {
                $ncut = 1;
            } else {
                $aa = substr($aa,1);
				$offset += 3;
            }
        }
    } elsif ($$ref{'stopcodon'} ne "yes" && $$ref{'startcodon'} eq "yes") {
        my $sub = '';
        my $ccut = 0; 
        while ($aa ne "" && $ccut == 0) {
            my $k = rindex($aa, 'K');
            my $r = rindex($aa, 'R');
            my $cut = $k > $r ? $k : $r;
            $cut = $cut > 0? $cut : -1; 
            my $sub = substr($aa, $cut+1);
            $aa = substr($aa,0,$cut+1);
            if ($aa eq "" || substr($sub,0,1) ne "P") {
                $ccut = 1;
            } else {
                $aa = substr($aa,0,$cut);
            }
        }
    } elsif ($$ref{'startcodon'} ne "yes" && $$ref{'stopcodon'} ne "yes") {
        my $sub = '';
        my $ccut = 0; 
        while ($aa ne "" && $ccut == 0) {
            my $k = rindex($aa, 'K');
            my $r = rindex($aa, 'R');
            my $cut = $k > $r ? $k : $r;
            $cut = $cut > 0? $cut : -1; 
            my $sub = substr($aa, $cut+1);
            $aa = substr($aa,0,$cut+1);
            if ($aa eq "" || substr($sub,0,1) ne "P") {
                $ccut = 1;
            } else {
                $aa = substr($aa,0,$cut);
            }
        }
        my $ncut = 0;
        while ($aa ne "" && $ncut == 0) {
            my $k = index($aa, 'K') >= 0? index($aa, 'K') : length($aa)-1;
            my $r = index($aa, 'R') >= 0? index($aa, 'R') : length($aa)-1;
            my $cut = $k < $r ? $k : $r;
            $aa = substr($aa, $cut+1);
	    $offset += 3*$cut+3;
            if ($aa eq "" || substr($aa,0,1) ne "P") {
                $ncut = 1;
            } else {
                $aa = substr($aa,1);
		$offset += 3;
            }
        }
    }
    my $currstart = 0;
    my $aacnt = 0;
    my @retaa = ();
    while ($aacnt < length($aa)){
		if (substr($aa,$aacnt,1) eq "*"){
			my $aac = substr($aa,$currstart,$aacnt+1-$currstart);
			$aac =~ tr/*/B/;
			if ($aac ne "B") {
				push @retaa,$aac;
				if ($$ref{'strand'} eq "+") {
					my $abspos = $$ref{'start'}+$offset+3*$currstart;
					my $refnuc = substr($seq,$offset+3*$currstart,1);
					print $o_tab qq($aac\t$gene\t$ctg\t$abspos\t$refnuc\t$refnuc\tassembly\n) if ($print);
				}else{
					my $abspos = $$ref{'stop'}-($offset+3*$currstart);
					my $refnuc = substr($seq,$offset+3*$currstart,1);
					$refnuc =~ tr/ACGT/TGCA/;
					print $o_tab qq($aac\t$gene\t$ctg\t$abspos\t$refnuc\t$refnuc\tassembly\n) if ($print);
				}
				last;
			}
		}elsif ((substr($aa,$aacnt,1) eq "K" || substr($aa,$aacnt,1) eq "R") && $aacnt+1 == length($aa)) {
			my $aac = substr($aa,$currstart,$aacnt+1-$currstart);
			push @retaa,$aac;
			if ($$ref{'strand'} eq "+") {
				my $abspos = $$ref{'start'}+$offset+3*$currstart;
				my $refnuc = substr($seq,$offset+3*$currstart,1);
				print $o_tab qq($aac\t$gene\t$ctg\t$abspos\t$refnuc\t$refnuc\tassembly\n) if ($print);
			}else{
				my $abspos = $$ref{'stop'}-($offset+3*$currstart);
				my $refnuc = substr($seq,$offset+3*$currstart,1);
				$refnuc =~ tr/ACGT/TGCA/;
				print $o_tab qq($aac\t$gene\t$ctg\t$abspos\t$refnuc\t$refnuc\tassembly\n) if ($print);
			}
		}elsif ( (substr($aa,$aacnt,1) eq "R" || substr($aa,$aacnt,1) eq "K") && substr($aa,$aacnt+1,1) ne "P"){
			my $aar = substr($aa,$currstart,$aacnt+1-$currstart);
			push @retaa,$aar;
			if ($$ref{'strand'} eq "+") {
				my $abspos = $$ref{'start'}+$offset+3*$currstart;
				my $refnuc = substr($seq,$offset+3*$currstart,1);
				print $o_tab qq($aar\t$gene\t$ctg\t$abspos\t$refnuc\t$refnuc\tassembly\n) if ($print);
			}else{
				my $abspos = $$ref{'stop'}-($offset+3*$currstart);
				my $refnuc = substr($seq,$offset+3*$currstart,1);
				$refnuc =~ tr/ACGT/TGCA/;
				print $o_tab qq($aar\t$gene\t$ctg\t$abspos\t$refnuc\t$refnuc\tassembly\n) if ($print);
			}
			$currstart = $aacnt +1;
		}elsif (substr($aa,$aacnt,1) ne "K" && substr($aa,$aacnt,1) ne "R" && $aacnt+1 == length($aa)){
			push @retaa,"";
		}
		$aacnt++;
    }
    unshift @retaa,$aa;
    return @retaa;
}

   

sub getAllVars {
    my ($ref, $variants) = @_;
    my @varcoll = ();
    my $vnum = 0;
    my %singleVar =();
    foreach my $v (sort {$$a{'pos'} <=> $$b{'pos'} || length($$b{'ref'}) <=> length($$a{'ref'})} @{$variants}){
		unless (exists $singleVar{$$v{'pos'}."-".$$v{'ref'}."-".$$v{'new'}}) {
			$singleVar{$$v{'pos'}."-".$$v{'ref'}."-".$$v{'new'}} = 1;
			push @varcoll, $v;
			$vnum++;
		}
    }
    my @newcoll = ();
    if ($vnum > 1) {
	my @testcol = [@varcoll];
	@newcoll = &getAdjVars(@testcol);
	my $i = 0;
	while ($newcoll[$i]){
	    my $ii = 0;
	    while ($newcoll[$ii]){
			my @var1 = ();
			my $j = 0;
			while ($newcoll[$i][$j]) {
				push @var1, $newcoll[$i][$j]{'pos'}.$newcoll[$i][$j]{'new'};
				$j++;
			}
			my $jj = 0;
			my @var2 = ();
			while ($newcoll[$ii][$jj]) {
				push @var2, $newcoll[$ii][$jj]{'pos'}.$newcoll[$ii][$jj]{'new'};
				$jj++;
			}
			my @minus = array_minus(@var1,@var2);
			if ($i != $ii && $#minus == -1) {
				splice @newcoll, $i, 1;
				$ii = 0;
			} else {
				$ii++;
			}
	    }
	    $i++;
	}
    } else {
	@newcoll = [@varcoll];
    }
    return @newcoll;
}


sub getAdjVars{
    my (@variants) = @_;
    my @varcoll = ();
    foreach my $v (@variants){
	push @varcoll,$v;
    }
    my @newcoll= ();
    my $i=0;
    while ($varcoll[$i]) {
	my $j = 0;
	while ($varcoll[$i][$j]) {
	    push @newcoll,$varcoll[$i][$j];
	    $j++;
	}
	$i++;
    }
    my @newcollA = [@newcoll];
    my $vid =0;
    while ($newcoll[$vid+1]){
		if ($newcoll[$vid]{'pos'}+length($newcoll[$vid]{'ref'})>$newcoll[$vid+1]{'pos'} ) {
			if ($newcoll[$vid]{'pos'}+length($newcoll[$vid]{'new'})>=$newcoll[$vid+1]{'pos'}+length($newcoll[$vid+1]{'new'}) && substr($newcoll[$vid]{'new'},$newcoll[$vid+1]{'pos'}-$newcoll[$vid]{'pos'},length($newcoll[$vid+1]{'new'})) eq $newcoll[$vid+1]{'new'}) {
				splice @newcoll, $vid+1, 1;
				shift @newcollA;
				push @newcollA, [@newcoll];
			}else{
				my @idx1 = (0..$#newcoll);
				my @idx2 = (0..$#newcoll);
				splice @idx1,$vid,1;
				splice @idx2,$vid+1,1;
				my @test1 = [@newcoll[@idx1]];
				my @newcoll1 = &getAdjVars(@test1);
				my @test2 = [@newcoll[@idx2]];
				my @newcoll2 = &getAdjVars(@test2);
				shift @newcollA;
				push @newcollA, @newcoll1;
				push @newcollA, @newcoll2;
				last;
			}
		} else {
			$vid++;
		}
    }
    return @newcollA;
}
 
  
sub getPeptides {
    my ($ctg,$ref, @variants) = @_;
    my $seq = "";
	my $offset = 0;
    my @vars = ();
    if ($$ref{'strand'} eq "+") {
		($seq,$offset) = &getSeqOnlyPlus($ctg,$$ref{'stop'} + 1 - $$ref{'start'},$$ref{'start'},$$ref{'startcodon'}); #get nucleotide sequence
		@vars = &reposVarPlus(@variants,$$ref{'start'}+$offset); #get position of variants relative to CDS
    } else{
		($seq,$offset) = &getSeqOnlyMinus($ctg,$$ref{'stop'}  + 1 - $$ref{'start'},$$ref{'start'},$$ref{'startcodon'});
		@vars = &reposVarMinus(@variants,$$ref{'stop'}-$offset);
    }
    #positions are now base 0!
    my %outpeps = ();
    my @peptideSeqs = &cutAndChange($seq,$$ref{'startcodon'},\@vars);
	if ($#peptideSeqs >= 0){
		if ($peptideSeqs[0] ne "too many") {
			foreach my $cseq (@peptideSeqs) {
				my @aapep = translatePeptide($$cseq{'seq'});
				foreach my $pep (@aapep){
					if ($pep eq "") {
						print $o_log qq(missing stop-codon in peptide of $$ref{'gene'}\n);
					} else{
						$outpeps{$pep} = $$cseq{'vars'};
					}
				}
			}
		} else{
			print $o_log qq(TOO MANY DIFFERENT VARIANTS IN PEPTIDES OF $$ref{'gene'}\n.);
		}
	}else{
		print $o_log qq(no cut peptides for a variant of $$ref{'gene'}\n);
	}
    return %outpeps;
}
 
sub getSeqOnlyPlus {
    my ($ctg,$length,$start,$startcodon) = @_;
    my $seq = substr($$seqbank{$ctg}, $start - 1, $length);
    $seq =~ tr/a-z/A-Z/;
    if ($length%3 !=0) {
		$seq = substr($seq,0,-($length%3));
    }
	my @codons = $seq =~ /(.{3})/g;
    my $aa = '';
    foreach my $cod (@codons) {
		my $add = $$aa_table{$cod} ? $$aa_table{$cod} : 'X';
		$aa .= $add;
		if ($add eq "*") {
			last
		}
    }
	my $offset = 0;
#cut peptides at K or R unless P
	if ($startcodon ne "yes") {
        my $ncut = 0;
        while ($aa ne "" && $ncut == 0) {
            my $k = index($aa, 'K') >= 0? index($aa, 'K') : length($aa)-1;
            my $r = index($aa, 'R') >= 0? index($aa, 'R') : length($aa)-1;
            my $cut = $k < $r ? $k : $r;
            $aa = substr($aa, $cut+1);
			$offset += 3*$cut+3;
			if ($aa eq "" || substr($aa,0,1) ne "P") {
                $ncut = 1;
            } else {
                $aa = substr($aa,1);
				$offset += 3;
            }
        }
    } 
    return ($seq,$offset);
} 

sub getSeqOnlyMinus {
    my ($ctg,$length,$start,$startcodon) = @_;
    my $seq = substr($$seqbank{$ctg}, $start - 1, $length);
    $seq =~ tr/a-z/A-Z/;
    $seq = join('', reverse(split(//, $seq)));
    $seq =~ tr/ACGT/TGCA/;
    if ($length%3 !=0) {
		$seq = substr($seq,0,-($length%3));
    }
		my @codons = $seq =~ /(.{3})/g;
    my $aa = '';
    foreach my $cod (@codons) {
		my $add = $$aa_table{$cod} ? $$aa_table{$cod} : 'X';
		$aa .= $add;
		if ($add eq "*") {
			last
		}
    }
	my $offset = 0;
#cut peptides at K or R unless P
	if ($startcodon ne "yes") {
        my $ncut = 0;
        while ($aa ne "" && $ncut == 0) {
            my $k = index($aa, 'K') >= 0? index($aa, 'K') : length($aa)-1;
            my $r = index($aa, 'R') >= 0? index($aa, 'R') : length($aa)-1;
            my $cut = $k < $r ? $k : $r;
            $aa = substr($aa, $cut+1);
			$offset += 3*$cut+3;
			if ($aa eq "" || substr($aa,0,1) ne "P") {
                $ncut = 1;
            } else {
                $aa = substr($aa,1);
				$offset += 3;
            }
        }
    } 
    return ($seq,$offset);
} 

sub reposVarPlus {
    my ($variants,$start) = @_;
    my @sortvar = ();
    foreach my $v (sort {$$b{'pos'} <=> $$a{'pos'}} @{$variants}){
		if($$v{'pos'}-$start >= 0){
			push @sortvar, {'pos' => $$v{'pos'}-$start, 'ref' => $$v{'ref'}, 'new' => $$v{'new'}, 'origpos' => $$v{'pos'}, 'origref' => $$v{'ref'}, 'orignew' => $$v{'new'}};
		}
	}	
    return @sortvar;
} 

sub reposVarMinus {
    my ($variants,$stop) = @_;
    my @sortvar = ();
    foreach my $v (sort {$$b{'pos'} <=> $$a{'pos'}} @{$variants}){
		if($stop-($$v{'pos'}+length($$v{'ref'})-1) >= 0){
			my $revref = join('', reverse(split(//, $$v{'ref'})));
			$revref =~ tr/ACGT/TGCA/;
			my $revnew = join('', reverse(split(//, $$v{'new'})));
			$revnew =~ tr/ACGT/TGCA/;
			push @sortvar, {'pos' => $stop-($$v{'pos'}+length($$v{'ref'})-1), 'ref' => $revref, 'new' => $revnew, 'origpos' => $$v{'pos'}, 'origref' => $$v{'ref'}, 'orignew' => $$v{'new'}};
		}
	}
    return @sortvar;
} 
   
sub cutAndChange{
    my ($seq,$start,$variants,$bgvariants) = @_;
    $seq =~ tr/a-z/A-Z/;
    my $lseq = $seq;
    my @codons = $seq =~ /(.{3})/g;
    my @variants = @{ $variants };
    my @bgvariants = ();
    my @sortvar = sort {$$b{'pos'} <=> $$a{'pos'}} @variants;
    #find if stopcodon was introduced in frameshift
    if ($start eq "fs") {
		my( @Sindex )= grep { $codons[$_] eq 'TAA' } 0..$#codons;
		push @Sindex,grep { $codons[$_] eq 'TAG' } 0..$#codons;
		push @Sindex,grep { $codons[$_] eq 'TGA' } 0..$#codons;
		if ($#Sindex >= 0) {
			my @sindex = sort { $a <=> $b } (@Sindex);
			foreach my $end (@sindex){
				my $endvar = 0;
				foreach my $v (@sortvar){
					if ($$v{'pos'} < 3*$end+3 && $$v{'pos'} + length($$v{'ref'})-1 >= 3*$end) {
						$endvar = 1;
					} 
				}
				if ($endvar == 0) {
					splice @codons, $end+1;
					$seq = substr($seq,0,3*$end+3);
					last;
				}
			}
		}
    }
    #find K and R codons
    my( @Rindex )= grep { $codons[$_] eq 'CGT' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'CGC' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'CGA' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'CGG' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'AGA' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'AGG' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'AAA' } 0..$#codons;
    push @Rindex,grep { $codons[$_] eq 'AAG' } 0..$#codons;
    #remove K and R cut positions followed by P
    my @cutindex = sort { $b <=> $a } (@Rindex);
    my $rind = $#cutindex;
    while ($rind >= 0) {
		if ($cutindex[$rind] != $#codons) {
			if($codons[$cutindex[$rind]+1] eq 'CCT' || $codons[$cutindex[$rind]+1] eq 'CCC' || $codons[$cutindex[$rind]+1] eq 'CCG'|| $codons[$cutindex[$rind]+1] eq 'CCA' ) {
				splice @cutindex,$rind,1;
			}
		}
		$rind--;
    }
    #sort cut positions (last first)
    push @cutindex, -1;
    #find/include variants per peptide
    my $lastpep = 1;
    my $oldcut = 0;
    my $subseq = "";
    my @keepseqs = ();
    my $firstpep = 0;
    my @allvars = ();
    my $cumvar = 0;
    my $startvar = $start eq "no" ? 0 : 1;
    foreach my $cut (@cutindex){
		if ($cut == -1) {
			$firstpep = 1;
		}
		my @currvar = ();
		my @cutvar = ();
		#get variants for this peptide
		foreach my $v (@sortvar){
			if ($lastpep == 1) {
				if($$v{'pos'} < length($seq) && $$v{'pos'} >= 3*$cut+3){
					push @currvar, $v;
				}elsif ($$v{'pos'} < 3*$cut+3 && $$v{'pos'} + length($$v{'ref'})-1 >= 3*$cut) {
					push @cutvar, $v;
				} 
				}elsif($firstpep == 1){
				if($$v{'pos'} < 3*$oldcut+3){
					push @currvar, $v;
				} 
				}else{
				if($$v{'pos'} < 3*$oldcut+3 && $$v{'pos'} >= 3*$cut+3){
					push @currvar, $v;
				}elsif ($$v{'pos'} < 3*$cut+3 && $$v{'pos'} + length($$v{'ref'})-1 >= 3*$cut) {
					push @cutvar, $v;
				} 
			}
		}
			#include variants
		my @currseqs = ();
		if ($lastpep == 1) {
			$subseq = substr($seq, 3*$cut+3);
		} else{
			$subseq = substr($seq, 3*$cut+3, 3*($oldcut-$cut));
		}
		my @evar = ();
		push @currseqs, {'seq' => $subseq, 'vars' => \@evar};
		while (@keepseqs){
			my $kept = shift @keepseqs;
			my $keptnew = $subseq.$$kept{'seq'};
			push @currseqs, {'seq' => $keptnew, 'vars' => $$kept{'vars'}};
		}
		if ($#currvar >= 0) {
			foreach my $v (@currvar) {
				if ($cumvar +  $#currseqs > $maxvar) {
					return "too many";
				}
				my $chvar = 0;
				if ($$v{'pos'} + length($$v{'ref'})-1 <= 2 && $firstpep == 1 && $startvar == 0){
					$startvar = 1;
					$chvar = 1;
				}
				if (length($$v{'ref'})==length($$v{'new'})) {
					my $seqnum = $#currseqs;
					my $seqcnt = 0;
					while ($seqcnt <= $seqnum){
						my $ss1 = substr($currseqs[$seqcnt]{'seq'}, 0, $$v{'pos'} - (3*$cut+3) );
						my $secstart2 = length($ss1) + length($$v{'ref'});
						my $ss2 ="";
						unless ($secstart2 >= length($currseqs[$seqcnt]{'seq'})) {
							$ss2 = substr($currseqs[$seqcnt]{'seq'}, $secstart2);
						}
						my $cvseq = $ss1.$$v{'new'}.$ss2;
						if (substr($cvseq, 0, 2) eq 'CC' && $firstpep == 0) {
							push @keepseqs, {'seq' => $cvseq, 'vars' => $currseqs[$seqcnt]{'vars'}};
							push @{$keepseqs[$#keepseqs]{'vars'}},$v;
						} else {
							push @currseqs, {'seq' => $cvseq, 'vars' => $currseqs[$seqcnt]{'vars'}};
							push @{$currseqs[$#keepseqs]{'vars'}},$v;
						}
						if ($firstpep == 1 && $startvar == 0) {
							my @tcodons = $cvseq =~ /(.{3})/g;
							my( @tRindex )= grep { $tcodons[$_] eq 'CGT' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'CGC' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'CGA' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'CGG' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'AGA' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'AGG' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'AAA' } 0..$#tcodons;
							push @tRindex,grep { $tcodons[$_] eq 'AAG' } 0..$#tcodons;
							#remove K and R cut positions followed by P
							my @tcutindex = sort { $b <=> $a } (@tRindex);
							my $trind = $#tcutindex;
							while ($trind >= 0) {
							if ($tcutindex[$trind] != $#tcodons) {
								if($tcodons[$tcutindex[$trind]+1] eq 'CCT' || $tcodons[$tcutindex[$trind]+1] eq 'CCC' || $tcodons[$tcutindex[$trind]+1] eq 'CCG'|| $tcodons[$tcutindex[$trind]+1] eq 'CCA' ) {
								splice @tcutindex,$trind,1;
								}
							}
							$trind--;
							}
							if ($#tcutindex>=0) {
								my $pseq = substr($cvseq,3*$tcutindex[$#tcutindex]+3);
								my @tseqs = (); 
								my @tevar = ();
								push @tseqs, {'seq' => $pseq, 'vars' => $currseqs[$seqcnt]{'vars'}};
								if($bgvariants) {
									push @{$tseqs[0]{'vars'}}, $bgvariants;
								}
								push @allvars, $tseqs[0];
							}
						}
						$seqcnt++;
					}
				}else{
					my $newSeq = substr($lseq, 3*$cut+3);
					my $ss1 = substr($newSeq, 0, $$v{'pos'} - (3*$cut+3) );
					my $secstart2 = length($ss1) + length($$v{'ref'});
					my $ss2 ="";
					unless ($secstart2 >= length($newSeq)) {
						$ss2 = substr($newSeq, $secstart2);
					}
					my $newseq = $ss1.$$v{'new'}.$ss2;
					if (length($newseq)%3 != 0) {
						$newseq = substr($newseq,0,-(length($newseq)%3));
					}
					my $altdiff = length($$v{'ref'})-length($$v{'new'});
					my @newVar = ();
					foreach my $subv (@sortvar){
						if ($$subv{'pos'} > $$v{'pos'}){
							push @newVar, {'pos' => $$subv{'pos'}-(3*$cut+3)-$altdiff, 'ref' => $$subv{'ref'}, 'new' => $$subv{'new'}, 'origpos' => $$subv{'origpos'}, 'origref' => $$subv{'origref'}, 'orignew' => $$subv{'orignew'}};
						}
					}
					
					if ($chvar == 1){
						$startvar = 0;
					}
					
					if (substr($newseq, 0, 2) eq 'CC' && $firstpep == 0) {
						my $currstart = "fs";
						my @peptideSeqs = &cutAndChange($newseq, $currstart,\@newVar,$v);
						push @keepseqs, @peptideSeqs;
					} elsif ($firstpep == 0){
						my $currstart = "fs";
						my @peptideSeqs = &cutAndChange($newseq, $currstart,\@newVar,$v);
						push @currseqs, @peptideSeqs;
					}else{
						my @peptideSeqs = &cutAndChange($newseq, $start,\@newVar,$v);
						push @currseqs, @peptideSeqs;
						$startvar = 1; 
					}
				}
			}
			
		}
		if ($#cutvar >= 0) {
			push @keepseqs, @currseqs;
		}
		$oldcut = $cut;
		if ($lastpep == 1) {
			$lastpep = 0;
		}
		unless ($startvar == 0 && $firstpep == 1) {
			foreach my $cseq (@currseqs){
				if($bgvariants) {
					push @{$$cseq{'vars'}}, $bgvariants;
				}
				push @allvars, $cseq;
			}
		}
		$cumvar += $#currseqs;
	}
    return @allvars;
} 

sub translatePeptide {
    my ($seq) = @_;
    $seq =~ tr/a-z/A-Z/;
    my @codons = $seq =~ /(.{3})/g;
    my $aa = '';
    foreach my $cod (@codons) {
	my $add = $$aab_table{$cod} ? $$aab_table{$cod} : 'X';
	$aa .= $add;
	if ($add eq "B") {
	    last
	}
    }
    my @retaa = ();
#cut peptides at K or R unless P
    my $currstart = 0;
    my $aacnt = 0;
    while ($aacnt < length($aa)){
	if (substr($aa,$aacnt,1) eq "B"){
	    my $aa = substr($aa,$currstart,$aacnt+1-$currstart);
	    push @retaa,$aa unless ($aa eq "B");
	    last;
	}elsif ((substr($aa,$aacnt,1) eq "K" || substr($aa,$aacnt,1) eq "R") && $aacnt+1 == length($aa)) {
	    push @retaa,substr($aa,$currstart,$aacnt+1-$currstart);
	}elsif ( (substr($aa,$aacnt,1) eq "K" || substr($aa,$aacnt,1) eq "R") && substr($aa,$aacnt+1,1) ne "P"){
	    push @retaa,substr($aa,$currstart,$aacnt+1-$currstart);
	    $currstart = $aacnt +1;
	}elsif (substr($aa,$aacnt,1) ne "K" && substr($aa,$aacnt,1) ne "R" && $aacnt+1 == length($aa)){
	    push @retaa,"";
	}
	$aacnt++;
    }
    return @retaa;
}

