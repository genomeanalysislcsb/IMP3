#!/usr/bin/perl

# script to rename all genes in a fasta file to MPA-suitable format (>generic|xxxxxxxx| genename)
# takes 1 input - the fasta file with the genes;
# writes new fasta file to standard out

# written by Anna Heintz-Buschart (November 2018, adapted March 2020)


use strict;
use warnings;
use Bio::DB::Fasta;

my $fastaFile = shift;

my $cnt = 1;
my $batch = "1";
my $geneNo = "";
my $db = Bio::DB::Fasta->new( $fastaFile );
my @ids = $db->get_all_ids;
foreach my $gene (@ids){
    my $sequence = $db->seq($gene);
    if  (!defined( $sequence )) {
            print STDERR "Sequence $gene not found. \n";
            next;
    }
    $geneNo = sprintf("%07d",$cnt);
    $cnt++;
    $gene = "generic|$batch$geneNo|$gene";
    print ">$gene\n", "$sequence\n";
    }   
}


sub get_all_ids {

 grep {!/^__/} keys %{shift->{offsets}}

}

