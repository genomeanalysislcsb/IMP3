#!/usr/bin/env perl

use Bio::DB::Fasta;

my $fastaFile = shift;
my $queryFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
open (IN, $queryFile);
while (<IN>){
    next if (/^\#/);
    chomp;
    my @f = split(/\t/);
    if ($f[2] eq "CDS") { #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
        my ($part) = ($f[8] =~ /partial=(\d\d);/);
        if($part eq "00"){
            my ($seq) = ($f[8] =~ /ID=([^;]+);/);
            my $sequence = $db->seq($seq);
            if  (!defined( $sequence )) {
                print STDERR "Sequence $seq not found. \n";
                next;
            }   
            print ">$seq\n", "$sequence\n";
        }
   }
}
close
