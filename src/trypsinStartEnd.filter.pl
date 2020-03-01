
#!/usr/bin/perl

use Bio::DB::Fasta;

my $fastaFile = $ARGV[0];
my $gffFile = $ARGV[1];
my $minPepLength = $ARGV[2];
my $minPepNum = $ARGV[3];

$minPepLength = $minPepLength ? $minPepLength : 6;
$minPepNum = $minPepNum ? $minPepNum : 2;

my $db = Bio::DB::Fasta->new( $fastaFile );


my $geneCount = 1;
open (IN, $gffFile);
while (<IN>) {
    next if (/^\#/);
    chomp;
    my @f = split(/\t/);
    if ($f[2] eq "CDS") { #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
        my $seq = ($f[8] =~ /ID=([^;]+);/);
        my $part = ($f[8] =~ /partial=(\d\d);/);
        if ($f[6] eq "-") {
            $part = reverse($part);
        }
        my $start = substr($part,0,1);
        my $end = substr($part,1,1);
        my $sequence = $db->seq($seq);
            if  (!defined( $sequence )) {
                print STDERR "Sequence $seq not found. \n";
                next;
            }
        if ($start ne "0" && $end eq "0") {
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
        } elsif ($end ne "0" && $start eq "0") {
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
        } elsif ($start ne "0" && $end ne "0") {
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
            $geneNo = sprintf("%07d",$geneCount);
            $geneCount++;
            my $batch = "1";
            $gene = "generic|$batch$geneNo|$seq";
            print ">$seq\n", "$sequence\n";
        }
    }
}
close(IN)
