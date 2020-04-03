#!/usr/bin/env perl

use strict;
use Bio::DB::Fasta;

my $fastaFile = shift;
my $gffFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );

my %loci =();
open (IN, $gffFile);
while (<IN>) {
    next if (/^\#/);
    chomp;
    my @f = split(/\t/);
				my $sequence = $db->seq($f[0],$f[3],$f[4]);
				defined( $sequence ) or die "Contig $seq not found in assembly. \n";
}
close(IN);
exit;

