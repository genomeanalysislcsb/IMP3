#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $infile = shift;

my $seq_io = Bio::SeqIO->new(-file => "$infile", -format => "fasta");

while (my $seq = $seq_io->next_seq){
    my $id = $seq->display_id;
    my $len = $seq->length;
    print "$id\t0\t$len\n";
    
}


