#!/usr/bin/env perl

use strict;
use Bio::DB::Fasta;

my $fastaFile = shift;
my $gffFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
my $aa_table = &loadAA();

my $gcnt = 1;
my $lcont = "";

my %loci =();
open (IN, $gffFile);
while (<IN>) {
    next if (/^\#/);
    chomp;
    my @f = split(/\t/);
    if ($f[2] eq "CDS") { #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
        my $seq = $f[0];
        if(exists($loci{$seq})){
         $loci{$seq} += 1;
        }else{
         $loci{$seq} = 1;
        }
        my $id = $seq."_".$loci{$seq};
        #print STDERR $seq;
        my ($part) = ($f[8] =~ /partial=(\d\d);/);
        if ($f[6] eq "-") {
            $part = reverse($part);
        }
        my $startC = substr($part,0,1);
	       my $strand = $f[6] eq "+" ? 1 : 2;
        my $sequence = $db->seq($seq,$f[3], $f[4], $strand);
            if  (!defined( $sequence )) {
                print STDERR "Sequence $seq not found. \n";
                next;
            }
	       my $pro = &translate($sequence);
	       substr($pro,0,1,"M") if $startC == 0;
        print qq(\>$id\n$pro\n);
    }
}
close(IN);
exit;

#### subroutines ###


sub translate {
    my $s = shift;
    my @codons = $s =~ /(.{3})/g;
    my $aa = '';
    foreach my $cod (@codons) {
	     my $add = $$aa_table{$cod} ? $$aa_table{$cod} : 'X';
	     $aa .= $add;
	     if ($add eq "*") {
	      last
	     }
    }
    return $aa;
}

sub loadAA {
    my %aa = ('ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'TTA' => 'L',
	      'TTG' => 'L', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TTT' => 'F', 'TTC' => 'F', 'ATG' => 'M',
	      'TGT' => 'C', 'TGC' => 'C', 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'GGT' => 'G', 'GGC' => 'G',
	      'GGA' => 'G', 'GGG' => 'G', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P', 'ACT' => 'T', 'ACC' => 'T',
	      'ACA' => 'T', 'ACG' => 'T', 'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	      'TAT' => 'Y', 'TAC' => 'Y', 'TGG' => 'W', 'CAA' => 'Q', 'CAG' => 'Q', 'AAT' => 'N', 'AAC' => 'N', 'CAT' => 'H',
	      'CAC' => 'H', 'GAA' => 'E', 'GAG' => 'E', 'GAT' => 'D', 'GAC' => 'D', 'AAA' => 'K', 'AAG' => 'K', 'CGT' => 'R',
	      'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*');
    return \%aa;
}


sub get_all_ids {

 grep {!/^__/} keys %{shift->{offsets}}

}
