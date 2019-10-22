#!/usr/bin/perl

# script to remove non-tryptic peptides from incomplete gene predictions
# takes 1 input: - the fasta file with the genes, as output by fixed MOCAT aux script;
# filters sequences to avoid empty or short sequences in the output
# writes new fasta file to standard out
# rule for tryptic peptides: "cut behind K or R unless followed by P"

# written by Anna Heintz-Buschart (January 2015 - bug fix in June 2016 - adapted to fasta-only input Oct 2017 - bug fix removed after the prodigal aux was correct, March 2018)
# there was a bug corrected in Feb 2018 in the Oct 2017 version, now it gives out the uncut complete sequences; latest version with filtering

use strict;
use Bio::DB::Fasta;

my $fastaFile = $ARGV[0];
my $minPepLength = $ARGV[1];
my $minPepNum = $ARGV[2];

$minPepLength = $minPepLength ? $minPepLength : 6;
$minPepNum = $minPepNum ? $minPepNum : 2;

my $db = Bio::DB::Fasta->new( $fastaFile );
my @ids = $db -> ids;

foreach my $id (@ids) {
    my $header = $db->header($id);
    my @f = split " ", $header;
    my $seq = $f[0];
    my $sequence = $db->seq($seq);
        if  (!defined( $sequence )) {
            print STDERR "Sequence $seq not found. \n";
            next;
        }
    (my $strand = $f[1]) =~ s/strand://g;
    (my $startc = $f[5]) =~ s/start_codon://g;
    (my $stopc = $f[6]) =~ s/stop_codon://g;
    if ($startc ne "yes" && $stopc eq "yes") {
        my $ncut = 0;
        while ($sequence ne "" && $ncut == 0) {
            my $k = index($sequence, 'K') >= 0? index($sequence, 'K') : length($sequence);
            my $r = index($sequence, 'R') >= 0? index($sequence, 'R') : length($sequence);
            my $cut = $k < $r ? $k : $r;
            $sequence = substr($sequence, $cut+1);
            if ($sequence eq "" || substr($sequence,0,1) ne "P") {
                $ncut = 1;
            } else {
                $sequence = substr($sequence,1);              
            }
        }
    } elsif ($stopc ne "yes" && $startc eq "yes") {
        my $sub = '';
        my $ccut = 0; 
        while ($sequence ne "" && $ccut == 0) {
            my $k = rindex($sequence, 'K');
            my $r = rindex($sequence, 'R');
            my $cut = $k > $r ? $k : $r;
            $cut = $cut > 0? $cut : -1; 
            my $sub = substr($sequence, $cut+1);
            $sequence = substr($sequence,0,$cut+1);
            if ($sequence eq "" || substr($sub,0,1) ne "P") {
                $ccut = 1;
            } else {
                $sequence = substr($sequence,0,$cut);
            }
        }
    } elsif ($startc ne "yes" && $stopc ne "yes") {
        my $sub = '';
        my $ccut = 0; 
        while ($sequence ne "" && $ccut == 0) {
            my $k = rindex($sequence, 'K');
            my $r = rindex($sequence, 'R');
            my $cut = $k > $r ? $k : $r;
            $cut = $cut > 0? $cut : -1; 
            my $sub = substr($sequence, $cut+1);
            $sequence = substr($sequence,0,$cut+1);
            if ($sequence eq "" || substr($sub,0,1) ne "P") {
                $ccut = 1;
            } else {
                $sequence = substr($sequence,0,$cut);
            }
        }
        my $ncut = 0;
        while ($sequence ne "" && $ncut == 0) {
            my $k = index($sequence, 'K') >= 0? index($sequence, 'K') : length($sequence);
            my $r = index($sequence, 'R') >= 0? index($sequence, 'R') : length($sequence);
            my $cut = $k < $r ? $k : $r;
            $sequence = substr($sequence, $cut+1);
            if ($sequence eq "" || substr($sequence,0,1) ne "P") {
                $ncut = 1;
            } else {
                $sequence = substr($sequence,1);              
            }
        }
    }
    
    my $currstart = 0;
    my $aacnt = 0;
    my @retaa = ();
    while ($aacnt < length($sequence)){
		if (substr($sequence,$aacnt,1) eq "*"){
			my $aac = substr($sequence,$currstart,$aacnt+1-$currstart);
			if ($aac ne "*") {
				push @retaa,$aac;
				last;
			}
		}elsif ((substr($sequence,$aacnt,1) eq "K" || substr($sequence,$aacnt,1) eq "R") && $aacnt+1 == length($sequence)) {
			my $aac = substr($sequence,$currstart,$aacnt+1-$currstart);
			push @retaa,$aac;
		}elsif ( (substr($sequence,$aacnt,1) eq "R" || substr($sequence,$aacnt,1) eq "K") && substr($sequence,$aacnt+1,1) ne "P"){
			my $aar = substr($sequence,$currstart,$aacnt+1-$currstart);
			push @retaa,$aar;
			$currstart = $aacnt +1;
		}elsif (substr($sequence,$aacnt,1) ne "K" && substr($sequence,$aacnt,1) ne "R" && $aacnt+1 == length($sequence)){
			push @retaa,"";
		}
		$aacnt++;
    }
    
    my $fpepcnt = 0;
	foreach my $tmppep (@retaa){
		if((length($tmppep)>$minPepLength && substr($tmppep,-1) eq "*") || (length($tmppep)>=$minPepLength && substr($tmppep,-1) ne "*")){
			$fpepcnt++;
		}
	}
    $sequence =~ s/\*//;
	if($fpepcnt >= $minPepNum){    
		print ">$seq\n", "$sequence\n";
	}
}

