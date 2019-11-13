#!/usr/bin/env perl

use strict;

my $file=$ARGV[0];

my %seqs=();
my @ids=();
my $id;

my $ct=0;
open(FILE,$file) or die $!;
while(my $str=<FILE>){
	chomp($str);
	$ct++;
	#next if length($str)==0;
	die "Empty string $ct\n" if length($str)==0;
	if ($str=~/^>(.+)$/){
		$id=$1;
		$seqs{$id}="";
	}else{
		$seqs{$id}.=$str;	
	}
}
close(FILE);

foreach $id (sort {$a<=>$b} keys %seqs){
    print $id."\t0\t".length($seqs{$id})."\n";
}

