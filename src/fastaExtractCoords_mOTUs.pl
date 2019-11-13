#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Fasta;
use Getopt::Long;

my ($fastaFile, $queryFile, $coordFile, $cogList, $outdir,$outname);
GetOptions('f|fasta=s' => \$fastaFile,
	   'm|motus=s' => \$queryFile,
       'l|location=s' => \$coordFile,
       'c|cogs=s' => \$cogList,
       'o|outdir=s' => \$outdir,
       'n|outname=s' => \$outname);

my $db = Bio::DB::Fasta->new( $fastaFile );

my @coga = split /,/, $cogList;
my %cogs = map { $_ => 1 } @coga;

my %seq_hash = ();

#meta_mOTU_v25_13437.metaMG0000001_COG0012       V1.FI14_revised_scaffold25112_1 101     1546
open (IN, $coordFile);
while (<IN>){
    chomp;
    my @fields = split /\t/; #0: motuName.ID_COG = complete name in fasta file 1: sequence name  2: first position  3: last position
    my $seq = $fields[0];
    my $spos = $fields[2];
    my $epos = $fields[3];
    my @att = split /\./, $seq;
    my $motu = $att[0];
    my @secatt = split /_/, $att[1];
    my $cog = $secatt[1];
    if(exists($cogs{$cog})){
        $seq_hash{$motu."_".$cog} = {'seq' => $seq, 'spos' => $spos, 'epos' => $epos};
    }
 }
close(IN);

my @outnames = map { $outdir."\/".$_.".".$outname.".fna" } @coga;

open (IN, $queryFile);
my %fhHash = map { open my $fh, '>', $_ ; $_ => $fh } @outnames;
while (<IN>){
    chomp;
    my @fields = split /\t/;
    my $motu = $fields[0];
    foreach my $cog (@coga){
        if(exists($seq_hash{$motu."_".$cog})){
            my $seq = $seq_hash{$motu."_".$cog}{'seq'};
            my $spos = $seq_hash{$motu."_".$cog}{'spos'};
            my $epos = $seq_hash{$motu."_".$cog}{'epos'};
            my $sequence = $db->seq($seq,$spos,$epos);
            next if (!defined( $sequence ));
            my $fhname = $outdir."\/".$cog.".".$outname.".fna" ;
           # print STDOUT $sequence;
            print {$fhHash{$fhname}} ">$seq\n$sequence\n";
    }}
}
close(IN);
map { close $_ } @outnames;
